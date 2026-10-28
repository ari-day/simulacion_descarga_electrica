
/*
 * Motor DBM 3D
 *
 * - Potencial eléctrico: ecuación de Laplace.
 * - Resolución numérica: SOR con estencil de 7 puntos.
 * - Crecimiento: selección probabilística de celdas candidatas.
 * - Campo eléctrico: E = -gradiente(V).
 * - Umbral atmosférico: Ec = E0 * delta * kHumedad.
 *
 * El modelo es una aproximación simplificada de la ramificación
 * de una descarga eléctrica. No reproduce todos los procesos
 * físicos de un rayo real.
 *
 * Referencias:
 * Niemeyer, Pietronero y Wiesmann (1984), PRL 52, 1033.
 * Sañudo et al. (1995), Nonlinear Processes in Geophysics, 2, 101-106.
 *
 * La corrección de humedad incluida es exploratoria y requiere
 * validación experimental o bibliográfica específica.
 */

(function (g) {
  "use strict";

  const R = 8.314462618;
  const MA = 0.0289647;
  const E0 = 3.0e6;
  const T0 = 293.15;
  const P0 = 101325;
  const U = 1.3;

  const LIBRE = 0;
  const CANAL = 1;
  const ELECTRODO = 2;

  const DEF = {
    nx: 33,
    ny: 33,
    nz: 33,
    eta: 1,
    tolerancia: 5e-4,
    maxIteraciones: 3000,
    maxCeldas: 1500,
    semilla: null,
    temperaturaC: 20,
    presionHpa: 1013.25,
    humedadRelativa: 50,
    dxMetros: 10,
    factorSobretension: U
  };

  const V6 = [
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1]
  ];

  // Generador pseudoaleatorio reproducible cuando se proporciona
  // una semilla. Sin semilla, utiliza Math.random().
  function aleatorio(s) {
    if (s == null) return Math.random;

    let a = s >>> 0;

    return () => {
      a = (a + 0x6D2B79F5) >>> 0;

      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);

      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Propiedades aproximadas del aire según las condiciones indicadas.
  function atmosfera(Tc, Ph, H) {
    const T = Tc + 273.15;
    const pPa = Ph * 100;

    const rho = pPa * MA / (R * T);
    const delta = (pPa / P0) * (T0 / T);

    // Presión de saturación aproximada, en hPa.
    const es = 6.1094 *
      Math.exp(17.625 * Tc / (Tc + 243.04));

    // Humedad absoluta aproximada, en g/m³.
    const ha = 216.7 * (H / 100) * es / T;

    // Corrección exploratoria de humedad.
    const hn = Math.max(
      1,
      Math.min(20, ha / Math.max(delta, 1e-9))
    );

    const kh = 1 + 0.010 * (hn - 11);

    return {
      densidad: rho,
      delta,
      kHumedad: kh,
      campoCritico: E0 * delta * kh
    };
  }

  function crear(o) {
    const p = Object.assign({}, DEF, o || {});

    const { nx, ny, nz } = p;

    if (![nx, ny, nz].every(
      v => Number.isInteger(v) && v >= 5
    )) {
      throw new Error(
        "nx, ny, nz deben ser enteros >= 5"
      );
    }

    if (!(p.eta >= 0)) {
      throw new Error("eta debe ser >= 0");
    }

    if (!(p.humedadRelativa >= 0 &&
          p.humedadRelativa <= 100)) {
      throw new Error("Humedad fuera de rango");
    }

    if (!(p.dxMetros > 0)) {
      throw new Error("dxMetros debe ser mayor que cero");
    }

    if (!(p.factorSobretension > 0)) {
      throw new Error(
        "factorSobretension debe ser mayor que cero"
      );
    }

    if (!(p.maxCeldas >= 1) ||
        !Number.isInteger(p.maxCeldas)) {
      throw new Error("maxCeldas debe ser un entero positivo");
    }

    if (!(p.maxIteraciones >= 1) ||
        !Number.isInteger(p.maxIteraciones)) {
      throw new Error(
        "maxIteraciones debe ser un entero positivo"
      );
    }

    if (!(p.tolerancia > 0)) {
      throw new Error("tolerancia debe ser mayor que cero");
    }

    const n = nx * ny * nz;
    const sl = nx * ny;

    const tipo = new Uint8Array(n);
    const phi = new Float64Array(n);

    /*
     * Potencial normalizado:
     * suelo: phi = 1
     * nube:  phi = 0
     *
     * El potencial físico se obtiene como Vfisico = phi * V.
     */

    for (let k = 0; k < sl; k++) {
      tipo[k] = ELECTRODO;
      phi[k] = 1;

      const superior = (nz - 1) * sl + k;

      tipo[superior] = ELECTRODO;
      phi[superior] = 0;
    }

    // Distribución inicial lineal del potencial.
    for (let z = 1; z < nz - 1; z++) {
      phi.fill(
        1 - z / (nz - 1),
        z * sl,
        (z + 1) * sl
      );
    }

    // Semilla inicial del canal en el centro del electrodo superior.
    const cx = nx >> 1;
    const cy = ny >> 1;
    const k0 = (nz - 1) * sl + cy * nx + cx;

    tipo[k0] = CANAL;
    phi[k0] = 0;

    const omega =
      2 / (1 + Math.sin(
        Math.PI / Math.max(nx, ny, nz)
      ));

    return {
      p,
      nx,
      ny,
      nz,
      tipo,
      phi,
      omega,
      rnd: aleatorio(p.semilla),
      estado: "listo",

      cel: [{
        x: cx,
        y: cy,
        z: nz - 1,
        k: k0,
        p: -1
      }],

      itU: 0,
      resU: 0,

      atm: atmosfera(
        p.temperaturaC,
        p.presionHpa,
        p.humedadRelativa
      ),

      // Tensión física de referencia entre los electrodos.
      V: p.factorSobretension *
        (nz - 1) * p.dxMetros * E0
    };
  }

  /*
   * Resuelve la ecuación de Laplace:
   *
   *     ∇²phi = 0
   *
   * Utiliza SOR y un estencil de 7 puntos.
   * Los bordes laterales se aproximan mediante reflexión,
   * correspondiente a una condición de Neumann.
   */
  function laplace(s) {
    const { nx, ny, nz, tipo, phi, omega } = s;
    const sl = nx * ny;

    const tol = s.p.tolerancia;
    const mx = s.p.maxIteraciones;

    let it = 0;
    let cm = Infinity;

    while (it < mx && cm > tol) {
      cm = 0;

      for (let z = 1; z < nz - 1; z++) {
        for (let y = 0; y < ny; y++) {
          const yb = y > 0 ? y - 1 : y + 1;
          const yf = y < ny - 1 ? y + 1 : y - 1;
          const o = z * sl + y * nx;

          for (let x = 0; x < nx; x++) {
            const k = o + x;

            if (tipo[k] !== LIBRE) continue;

            const xl = x > 0 ? x - 1 : x + 1;
            const xr = x < nx - 1 ? x + 1 : x - 1;

            const pr = (
              phi[o + xl] +
              phi[o + xr] +
              phi[z * sl + yb * nx + x] +
              phi[z * sl + yf * nx + x] +
              phi[k + sl] +
              phi[k - sl]
            ) / 6;

            const d = omega * (pr - phi[k]);

            phi[k] += d;

            const ad = Math.abs(d);

            if (ad > cm) cm = ad;
          }
        }
      }

      it++;
    }

    s.itU = it;
    s.resU = cm;

    return cm <= tol;
  }

  /*
   * Calcula el campo eléctrico local a partir del gradiente
   * del potencial normalizado.
   *
   * E = -gradiente(Vfisico)
   * Vfisico = phi * s.V
   *
   * En los bordes laterales se impone componente normal nula,
   * de acuerdo con la aproximación de Neumann utilizada
   * por el solucionador de Laplace.
   *
   * Se supone una separación uniforme entre celdas en los
   * tres ejes: dxMetros.
   */
  function campoLocal(s, k) {
    const { nx, ny, nz, phi, p } = s;
    const sl = nx * ny;

    const x = k % nx;
    const y = Math.floor(k / nx) % ny;
    const z = Math.floor(k / sl);

    let dphidx = 0;
    let dphidy = 0;
    let dphidz = 0;

    // Derivada en X: diferencia centrada en el interior.
    // En el borde lateral, la condición de Neumann
    // implica gradiente normal aproximadamente nulo.
    if (x > 0 && x < nx - 1) {
      dphidx = (phi[k + 1] - phi[k - 1]) / 2;
    }

    // Derivada en Y.
    if (y > 0 && y < ny - 1) {
      dphidy = (phi[k + nx] - phi[k - nx]) / 2;
    }

    // Derivada en Z.
    // Las candidatas de crecimiento están entre los electrodos.
    if (z > 0 && z < nz - 1) {
      dphidz = (phi[k + sl] - phi[k - sl]) / 2;
    } else if (z === 0) {
      dphidz = phi[k + sl] - phi[k];
    } else {
      dphidz = phi[k] - phi[k - sl];
    }

    // Convierte el gradiente normalizado a V/m.
    const escala = s.V / p.dxMetros;

    const Ex = -escala * dphidx;
    const Ey = -escala * dphidy;
    const Ez = -escala * dphidz;

    return Math.sqrt(
      Ex * Ex + Ey * Ey + Ez * Ez
    );
  }

  /*
   * Incorpora una nueva celda al canal mediante una
   * selección probabilística basada en el campo local.
   */
  function paso(s) {
    if (s.estado !== "listo" &&
        s.estado !== "creciendo") {
      return null;
    }

    const {
      nx, ny, nz, tipo, phi, cel, p
    } = s;

    if (cel.length >= p.maxCeldas) {
      s.estado = "limite";
      return null;
    }

    if (s.estado === "listo") {
      if (!laplace(s)) {
        s.estado = "error_convergencia";
        return null;
      }

      s.estado = "creciendo";
    }

    const vis = new Set();
    const cand = [];
    const w = [];

    const Ec = s.atm.campoCritico;

    let suma = 0;

    // Examina las seis vecinas de cada celda del canal.
    for (let i = 0; i < cel.length; i++) {
      for (const [a, b, d] of V6) {
        const x = cel[i].x + a;
        const y = cel[i].y + b;
        const z = cel[i].z + d;

        if (
          x < 0 || y < 0 || z < 0 ||
          x >= nx || y >= ny || z >= nz
        ) {
          continue;
        }

        const k = (z * ny + y) * nx + x;

        if (tipo[k] !== LIBRE || vis.has(k)) {
          continue;
        }

        vis.add(k);

        // CORRECCIÓN:
        // ahora se calcula el campo mediante el gradiente
        // del potencial, no con abs(phi[k]) * V/dx.
        const Ef = campoLocal(s, k);

        // Regla de crecimiento con umbral crítico.
        if (Ef > Ec) {
          const wi = Math.pow(
            Ef / Ec - 1,
            p.eta
          );

          if (wi > 0 && Number.isFinite(wi)) {
            cand.push({
              x,
              y,
              z,
              k,
              p: i
            });

            w.push(wi);
            suma += wi;
          }
        }
      }
    }

    if (!cand.length || !(suma > 0)) {
      s.estado = cel.length === 1
        ? "ruptura_no_iniciada"
        : "sin_candidatas";

      return null;
    }

    // Selección aleatoria ponderada.
    let r = s.rnd() * suma;
    let j = 0;

    for (; j < w.length - 1; j++) {
      r -= w[j];

      if (r < 0) break;
    }

    const e = cand[j];

    tipo[e.k] = CANAL;
    phi[e.k] = 0;

    cel.push(e);

    // Recalcula el potencial después de modificar el canal.
    if (!laplace(s)) {
      s.estado = "error_convergencia";
      return e;
    }

    if (e.z === 1) {
      s.estado = "electrodo_alcanzado";
    } else if (cel.length >= p.maxCeldas) {
      s.estado = "limite";
    }

    return e;
  }

  // Estimación complementaria de dimensión fractal por conteo de cajas.
  function dimensionFractal(s) {
    if (s.cel.length < 20) return null;

    const xs = [];
    const ys = [];

    for (
      let e = 1;
      e <= Math.min(s.nx, s.ny, s.nz) / 4;
      e *= 2
    ) {
      const c = new Set();

      for (const q of s.cel) {
        c.add(
          Math.floor(q.x / e) + "," +
          Math.floor(q.y / e) + "," +
          Math.floor(q.z / e)
        );
      }

      xs.push(Math.log(1 / e));
      ys.push(Math.log(c.size));
    }

    if (xs.length < 3) return null;

    const m = a =>
      a.reduce((u, v) => u + v, 0) / a.length;

    const mx = m(xs);
    const my = m(ys);

    let sxy = 0;
    let sxx = 0;

    xs.forEach((x, i) => {
      sxy += (x - mx) * (ys[i] - my);
      sxx += (x - mx) ** 2;
    });

    return sxx > 0 ? sxy / sxx : null;
  }

  const api = {
    constantes: {
      E0,
      U,
      LIBRE,
      CANAL,
      ELECTRODO
    },
    crear,
    paso,
    laplace,
    campoLocal,
    atmosfera,
    dimensionFractal
  };

  g.MotorDBM3D = api;

  if (
    typeof module !== "undefined" &&
    module.exports
  ) {
    module.exports = api;
  }
})(typeof window !== "undefined" ? window : globalThis);
