// createStreamline.js
(function (window) {
  /**
   * Main entry point.
   * @param {Array<number>} x - 1D, evenly spaced x-domain array.
   * @param {Array<number>} y - 1D, evenly spaced y-domain array.
   * @param {Array<Array<number>>} u - 2D array for the x-component of the vector field.
   * @param {Array<Array<number>>} v - 2D array for the y-component of the vector field.
   * @param {Object} options - Optional parameters.
   *        options.density     - Controls streamline density (default 1).
   *        options.angle       - Arrowhead angle in radians (default Math.PI/9).
   *        options.arrow_scale - Arrowhead scale factor (default 0.09).
   *        options.maxLength   - Maximum integration length (default 2).
   * @returns {Object} An object with properties { x, y } containing the streamline and arrow data.
   */
  function createStreamline(x, y, u, v, options) {
    options = options || {};
    const density = options.density !== undefined ? options.density : 1;
    const angle = options.angle !== undefined ? options.angle : Math.PI / 9;
    const arrow_scale =
      options.arrow_scale !== undefined ? options.arrow_scale : 0.09;
    const maxLength =
      options.maxLength !== undefined ? options.maxLength : 2; // new parameter

    // (Optional) Validate that x and y are evenly spaced.
    // For brevity, we assume they are.

    const stream = new _Streamline(x, y, u, v, density, angle, arrow_scale, maxLength);
    // Sum all streamline trajectories
    const streamlineData = stream.sumStreamlines();
    // Get arrow head coordinates
    const arrowData = stream.getStreamlineArrows();
    // Concatenate streamline and arrow data
    const finalX = streamlineData.x.concat(arrowData.x);
    const finalY = streamlineData.y.concat(arrowData.y);
    return { x: finalX, y: finalY };
  }

  // -----------------------------------------------------
  // _Streamline class – handles integration and arrow creation
  // -----------------------------------------------------
  class _Streamline {
    constructor(x, y, u, v, density, angle, arrow_scale, maxLength) {
      this.x = x; // 1D array (domain)
      this.y = y; // 1D array (domain)
      this.u = u; // 2D array, dimensions: [y.length][x.length]
      this.v = v; // 2D array, same dimensions as u
      this.angle = angle;
      this.arrow_scale = arrow_scale;
      this.density = Math.floor(30 * density); // Scale density similarly to Python.
      this.delta_x = x[1] - x[0];
      this.delta_y = y[1] - y[0];
      this.maxLength = maxLength; // new max integration length parameter

      // Set up spacing for seed grid:
      this.spacing_x = x.length / (this.density - 1);
      this.spacing_y = y.length / (this.density - 1);
      // Create blank grid to track used seeds.
      this.blank = Array.from({ length: this.density }, () =>
        Array(this.density).fill(0)
      );
      this.trajectories = [];

      // Rescale u and v onto index space.
      const domainWidth = x[x.length - 1] - x[0];
      const domainHeight = y[y.length - 1] - y[0];
      this.u_scaled = [];
      this.v_scaled = [];
      this.speed = [];
      for (let j = 0; j < this.v.length; j++) {
        this.u_scaled[j] = [];
        this.v_scaled[j] = [];
        this.speed[j] = [];
        for (let i = 0; i < this.u[j].length; i++) {
          const u_val = (this.u[j][i] / domainWidth) * x.length;
          const v_val = (this.v[j][i] / domainHeight) * y.length;
          this.u_scaled[j][i] = u_val;
          this.v_scaled[j][i] = v_val;
          this.speed[j][i] = Math.sqrt(u_val * u_val + v_val * v_val);
        }
      }

      this.st_x = [];
      this.st_y = [];
      this.getStreamlines();
    }

    // Given continuous positions xi, yi (in index space), return grid indices.
    blankPos(xi, yi) {
      const pos_x = Math.floor(xi / this.spacing_x + 0.5);
      const pos_y = Math.floor(yi / this.spacing_y + 0.5);
      return [pos_x, pos_y];
    }

    // Bilinear interpolation on a 2D array 'a' at continuous indices (xi, yi)
    valueAt(a, xi, yi) {
      const x0 = Math.floor(xi);
      const y0 = Math.floor(yi);
      const x1 = x0 + 1;
      const y1 = y0 + 1;
      // Safeguard: if indices are out-of-bound, return 0.
      if (
        x0 < 0 ||
        y0 < 0 ||
        x1 >= this.x.length ||
        y1 >= this.y.length ||
        !a[y0] ||
        !a[y1]
      ) {
        return 0;
      }
      const a00 = a[y0][x0],
        a01 = a[y0][x1],
        a10 = a[y1][x0],
        a11 = a[y1][x1];
      const xt = xi - x0;
      const yt = yi - y0;
      const a0 = a00 * (1 - xt) + a01 * xt;
      const a1 = a10 * (1 - xt) + a11 * xt;
      return a0 * (1 - yt) + a1 * yt;
    }

    // Check if continuous indices xi, yi are within valid bounds.
    check(xi, yi) {
      return (
        xi >= 0 &&
        xi < this.x.length - 1 &&
        yi >= 0 &&
        yi < this.y.length - 1
      );
    }

    // RK4 integration with added checks and using maxLength.
    rk4(x0, y0, f) {
      const ds = 0.01;
      let s_total = 0;
      let xi = x0;
      let yi = y0;
      const [init_bx, init_by] = this.blankPos(xi, yi);
      let current_bx = init_bx,
        current_by = init_by;
      const traj_x = [xi];
      const traj_y = [yi];
      const blank_changes = [];
      while (this.check(xi, yi)) {
        const k1 = f(xi, yi);
        const mid1_x = xi + 0.5 * ds * k1[0];
        const mid1_y = yi + 0.5 * ds * k1[1];
        if (!this.check(mid1_x, mid1_y)) break;
        const k2 = f(mid1_x, mid1_y);
        const mid2_x = xi + 0.5 * ds * k2[0];
        const mid2_y = yi + 0.5 * ds * k2[1];
        if (!this.check(mid2_x, mid2_y)) break;
        const k3 = f(mid2_x, mid2_y);
        const mid3_x = xi + ds * k3[0];
        const mid3_y = yi + ds * k3[1];
        if (!this.check(mid3_x, mid3_y)) break;
        const k4 = f(mid3_x, mid3_y);
        const delta_x =
          (ds * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])) / 6.0;
        const delta_y =
          (ds * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])) / 6.0;
        xi += delta_x;
        yi += delta_y;
        if (!this.check(xi, yi)) break;
        s_total += ds;
        // Use the new maxLength parameter:
        if (s_total > this.maxLength) break;
        traj_x.push(xi);
        traj_y.push(yi);
        const [new_bx, new_by] = this.blankPos(xi, yi);
        if (new_bx !== current_bx || new_by !== current_by) {
          if (this.blank[new_by] && this.blank[new_by][new_bx] === 0) {
            this.blank[new_by][new_bx] = 1;
            blank_changes.push([new_bx, new_by]);
            current_bx = new_bx;
            current_by = new_by;
          } else {
            break;
          }
        }
      }
      return { s_total, traj_x, traj_y, blank_changes };
    }

    // Integrate a trajectory in both forward and backward directions.
    rk4Integrate(x0, y0) {
      const f = (xi, yi) => {
        const speed_val = this.valueAt(this.speed, xi, yi);
        const dt_ds = speed_val !== 0 ? 1.0 / speed_val : 0;
        const ui = this.valueAt(this.u_scaled, xi, yi);
        const vi = this.valueAt(this.v_scaled, xi, yi);
        return [ui * dt_ds, vi * dt_ds];
      };
      const g = (xi, yi) => {
        const speed_val = this.valueAt(this.speed, xi, yi);
        const dt_ds = speed_val !== 0 ? 1.0 / speed_val : 0;
        const ui = this.valueAt(this.u_scaled, xi, yi);
        const vi = this.valueAt(this.v_scaled, xi, yi);
        return [-ui * dt_ds, -vi * dt_ds];
      };
      const forward = this.rk4(x0, y0, f);
      const backward = this.rk4(x0, y0, g);
      const s_total = forward.s_total + backward.s_total;
      const traj_x = backward.traj_x.reverse().concat(forward.traj_x.slice(1));
      const traj_y = backward.traj_y.reverse().concat(forward.traj_y.slice(1));
      if (traj_x.length < 1) return null;
      if (s_total > 0.2) {
        const [init_bx, init_by] = this.blankPos(x0, y0);
        this.blank[init_by][init_bx] = 1;
        return { traj_x, traj_y };
      } else {
        return null;
      }
    }

    // Attempt to compute a trajectory from a given seed.
    traj(xb, yb) {
      if (xb < 0 || xb >= this.density || yb < 0 || yb >= this.density) return;
      if (this.blank[yb][xb] === 0) {
        const x0 = xb * this.spacing_x;
        const y0 = yb * this.spacing_y;
        const result = this.rk4Integrate(x0, y0);
        if (result !== null) {
          this.trajectories.push(result);
        }
      }
    }

    // Compute streamlines by scanning a grid with indentations.
    getStreamlines() {
      for (let indent = 0; indent < Math.floor(this.density / 2); indent++) {
        for (let xi = 0; xi < this.density - 2 * indent; xi++) {
          this.traj(xi + indent, indent);
          this.traj(xi + indent, this.density - 1 - indent);
          this.traj(indent, xi + indent);
          this.traj(this.density - 1 - indent, xi + indent);
        }
      }
      for (let t = 0; t < this.trajectories.length; t++) {
        const traj = this.trajectories[t];
        const conv_x = traj.traj_x.map(val => val * this.delta_x + this.x[0]);
        const conv_y = traj.traj_y.map(val => val * this.delta_y + this.y[0]);
        conv_x.push(NaN);
        conv_y.push(NaN);
        this.st_x = this.st_x || [];
        this.st_y = this.st_y || [];
        this.st_x.push(conv_x);
        this.st_y.push(conv_y);
      }
    }

    // Create arrow head coordinates for each streamline.
    getStreamlineArrows() {
      const arrows_x = [];
      const arrows_y = [];
      for (let i = 0; i < this.st_x.length; i++) {
        const xs = this.st_x[i];
        const ys = this.st_y[i];
        if (xs.length < 3) continue;
        const idx = Math.floor(xs.length / 3);
        const arrow_end_x = xs[idx];
        const arrow_start_x = xs[idx - 1];
        const arrow_end_y = ys[idx];
        const arrow_start_y = ys[idx - 1];
        const dif_x = arrow_end_x - arrow_start_x;
        const dif_y = arrow_end_y - arrow_start_y;
        const streamline_ang = Math.atan2(dif_y, dif_x);
        const ang1 = streamline_ang + this.angle;
        const ang2 = streamline_ang - this.angle;
        const seg1_x = Math.cos(ang1) * this.arrow_scale;
        const seg1_y = Math.sin(ang1) * this.arrow_scale;
        const seg2_x = Math.cos(ang2) * this.arrow_scale;
        const seg2_y = Math.sin(ang2) * this.arrow_scale;
        let point1_x, point1_y, point2_x, point2_y;
        if (dif_x >= 0) {
          point1_x = arrow_end_x - seg1_x;
          point1_y = arrow_end_y - seg1_y;
          point2_x = arrow_end_x - seg2_x;
          point2_y = arrow_end_y - seg2_y;
        } else {
          point1_x = arrow_end_x + seg1_x;
          point1_y = arrow_end_y + seg1_y;
          point2_x = arrow_end_x + seg2_x;
          point2_y = arrow_end_y + seg2_y;
        }
        arrows_x.push(point1_x, arrow_end_x, point2_x, NaN);
        arrows_y.push(point1_y, arrow_end_y, point2_y, NaN);
      }
      return { x: arrows_x, y: arrows_y };
    }

    // Concatenate all streamlines into a single trace.
    sumStreamlines() {
      let all_x = [];
      let all_y = [];
      for (let i = 0; i < this.st_x.length; i++) {
        all_x = all_x.concat(this.st_x[i]);
        all_y = all_y.concat(this.st_y[i]);
      }
      return { x: all_x, y: all_y };
    }
  }

  window.createStreamline = createStreamline;
})(window);
