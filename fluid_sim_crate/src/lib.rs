#![no_std]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

pub mod FluidSimulation {

    use libm::{floorf, sqrtf};

    static max_particles_setting: usize = 400;
    static number_of_vertical_cells_setting: usize = 23;
    static number_of_horizontal_cells_setting: usize = 23;
    static max_particles_x2_setting: usize = max_particles_setting * 2;
    static number_of_cells_setting: usize =
        number_of_vertical_cells_setting * number_of_horizontal_cells_setting;
    static number_of_cells_x2_setting: usize = number_of_cells_setting * 2;
    static number_of_cells_setting_plus1: usize = number_of_cells_setting + 1;

    static simHeight: f32 = 23.0;
    static simWidth: f32 = 23.0;

    #[derive(PartialEq, Clone, Copy)]
    pub enum CellType {
        AIR_CELL,
        FLUID_CELL,
        SOLID_CELL,
    }

    pub struct FlipFluid {
        density: f32,

        /// the number of cells in the x direction
        fNumX: f32,

        /// the number of cells in the y direction
        fNumY: f32,

        /// the largest distance between 2 adjacent cells (basically the size of a cell)
        h: f32,

        /// the inverse of the largest distance between 2 adjacent cells (1.0/h)
        fInvSpacing: f32,

        /// the number of total cells
        fNumCells: f32,

        /// an array of the cell positions 2*i is y coordinate, 2*i+1 is x coordinate
        u: [f32; number_of_cells_x2_setting],

        /// an array of the cell velocities 2*i is the vertical velocity, 2*i+1 is the horizontal velocity
        v: [f32; number_of_cells_x2_setting],

        /// this is an array of the amount the position will change in one time step 2*i is the
        /// y coordinate change and 2*i+1 is the x coordinate change
        du: [f32; number_of_cells_x2_setting],

        /// this is an array of the amount the velocity will change in one time step 2*i is the
        /// vertical velocity change and 2*i+1 is the horizontal velocity change.
        dv: [f32; number_of_cells_x2_setting],

        /// this is an array of the previous iteration's particle position 2*i are y coordinates
        /// and 2*i+1 are x coordinates
        prevU: [f32; number_of_cells_x2_setting],

        /// this is an array of the previous iteration's particle velocity.  2*i are y vertical
        /// velocities and 2*i+1 are horizontal velocities.
        prevV: [f32; number_of_cells_x2_setting],

        p: [f32; number_of_cells_setting],
        s: [f32; number_of_cells_setting],
        cellType: [CellType; number_of_cells_setting],

        /// the max number of particles, used to size the arrays
        _maxParticles: i32,

        /// the particle position, in meters, indexes 2*i are vertical, 2*i+1 are horizontal
        pub particlePos: [f32; max_particles_x2_setting],

        /// the particle velocity, in m/s, indexes 2*i are vertical, 2*i+1 are horizontal
        particleVel: [f32; max_particles_x2_setting],

        particleDensity: [f32; number_of_cells_setting],
        particleRestDensity: f32,
        particleRadius: f32,

        /// pInvSpacing = 1.0 / (2.2 * particleRadius);
        pInvSpacing: f32,

        /// the number of cells in the x direction
        pNumX: i32,

        /// the number of cells in the y direction
        pNumY: i32,

        /// the total number of cells in the array
        pNumCells: i32,

        numCellParticles: [i32; number_of_cells_setting],
        firstCellParticle: [i32; number_of_cells_setting_plus1],
        cellParticleIds: [i32; max_particles_setting],

        /// the total number of particles in the system
        pub numParticles: i32,
    }

    impl FlipFluid {
        fn new(
            density: f32,
            width: f32,
            height: f32,
            spacing: f32,
            particleRadius: f32,
            _maxParticles: i32,
        ) -> FlipFluid {
            // this.density = density;
            // density
            // this.fNumX = Math.floor(width / spacing) + 1;
            let fNumX = floorf(width / spacing);
            // this.fNumY = Math.floor(height / spacing) + 1;
            let fNumY = floorf(height / spacing);
            // this.h = Math.max(width / this.fNumX, height / this.fNumY);
            let h = (width / fNumX).max(height / fNumY);
            // this.fInvSpacing = 1.0 / this.h;
            let fInvSpacing = 1.0 / h;
            // this.fNumCells = this.fNumX * this.fNumY;
            let fNumCells = fNumX * fNumY;

            // this.u = new Float32Array(this.fNumCells);
            let u = [0f32; number_of_cells_x2_setting];
            // this.v = new Float32Array(this.fNumCells);
            let v = [0f32; number_of_cells_x2_setting];
            // this.du = new Float32Array(this.fNumCells);
            let du = [0f32; number_of_cells_x2_setting];
            // this.dv = new Float32Array(this.fNumCells);
            let dv = [0f32; number_of_cells_x2_setting];
            // this.prevU = new Float32Array(this.fNumCells);
            let prevU = [0f32; number_of_cells_x2_setting];
            // this.prevV = new Float32Array(this.fNumCells);
            let prevV = [0f32; number_of_cells_x2_setting];
            // this.p = new Float32Array(this.fNumCells);
            let p = [0f32; number_of_cells_setting];
            // this.s = new Float32Array(this.fNumCells);
            let s = [0f32; number_of_cells_setting];
            // this.cellType = new Int32Array(this.fNumCells);
            let cellType = [CellType::AIR_CELL; number_of_cells_setting];
            // this.particleRadius = particleRadius;
            // this.maxParticles = maxParticles;
            //maxParticles
            // this.particlePos = new Float32Array(2 * this.maxParticles);
            let mut particlePos = [0.0; max_particles_x2_setting];
            let mut count: usize = 0;
            for i in 1..21 {
                for j in 1..21 {
                    particlePos[count * 2] = j as f32;
                    particlePos[count * 2 + 1] = i as f32;
                    count += 1;
                }
            }
            //

            // this.particleVel = new Float32Array(2 * this.maxParticles);
            let particleVel = [0.0; max_particles_x2_setting];
            // this.particleDensity = new Float32Array(this.fNumCells);
            let particleDensity = [0f32; number_of_cells_setting];
            // this.particleRestDensity = 0.0;
            let particleRestDensity = 0.0f32;

            //particleRadius
            // this.pInvSpacing = 1.0 / (2.2 * particleRadius);
            // let pInvSpacing = 1.0 / (2.2 * particleRadius); // original
            let pInvSpacing = 1.0;
            // this.pNumX = Math.floor(width * this.pInvSpacing) + 1;
            let pNumX = floorf(width * pInvSpacing) as i32;
            // this.pNumY = Math.floor(height * this.pInvSpacing) + 1;
            let pNumY = floorf(height * pInvSpacing) as i32;
            // this.pNumCells = this.pNumX * this.pNumY;
            let pNumCells = pNumX * pNumY;
            // this.numCellParticles = new Int32Array(this.pNumCells);
            let numCellParticles = [0; number_of_cells_setting];
            // this.firstCellParticle = new Int32Array(this.pNumCells + 1);
            let firstCellParticle = [0; number_of_cells_setting_plus1];
            // this.cellParticleIds = new Int32Array(maxParticles);
            let cellParticleIds = [0; max_particles_setting];
            // this.numParticles = 0;
            let numParticles = 0i32;

            FlipFluid {
                density,
                fNumX,
                fNumY,
                h,
                fInvSpacing,
                fNumCells,
                u,
                v,
                du,
                dv,
                prevU,
                prevV,
                p,
                s,
                cellType,
                _maxParticles,
                particlePos,
                particleVel,
                particleDensity,
                particleRestDensity,
                particleRadius,
                pInvSpacing,
                pNumY,
                numCellParticles,
                cellParticleIds,
                pNumX,
                pNumCells,
                firstCellParticle,
                numParticles,
            }
        }

        // integrateParticles(dt, gravity)
        /// add gravity to the velocity of the particles and calculate positions.
        ///
        /// for velocity (vertical only): Vert_Vel_new = Vert_Vel_old + dt * acceleration
        ///
        /// vertical velocity is unchanged here Horiz_Vel_New = Horiz_Vel_Old
        ///
        /// the position of each is then calculated
        ///
        /// Vert_Pos_new = Vert_Pos_old + Vert_Vel_New * dt
        ///
        /// Horiz_Pos_new = Horz_Pos_old + Horiz_Vel_New * dt
        fn integrateParticles(&mut self, dt: f32, yGravity: f32, xGravity: f32) {
            // for (var i = 0; i < this.numParticles; i++) {

            for i in 0..self.numParticles {
                self.particleVel[(2 * i) as usize] += dt * xGravity;
                // this.particleVel[2 * i + 1] += dt * gravity;
                self.particleVel[(2 * i + 1) as usize] += dt * yGravity;
                // this.particlePos[2 * i] += this.particleVel[2 * i] * dt;
                self.particlePos[(2 * i) as usize] += self.particleVel[(2 * i) as usize] * dt;
                // this.particlePos[2 * i + 1] += this.particleVel[2 * i + 1] * dt;
                self.particlePos[(2 * i + 1) as usize] +=
                    self.particleVel[(2 * i + 1) as usize] * dt;
            }
        }

        fn showParticles(&mut self) {
            for i in 0..self.numParticles {
                let mut cell_location: (usize, usize) = (0, 0);
                cell_location.0 = floorf(self.particlePos[2 * i as usize]) as usize;
                cell_location.1 = floorf(self.particlePos[2 * i as usize + 1]) as usize;
                self.cellType[cell_location.1 * 23 + cell_location.0] = CellType::FLUID_CELL;
            }
        }
        // pushParticlesApart(numIters)
        /// store the particle positions in x and y and
        /// make incompressible by making sure the amount of fluid that enters a cell is equal to the fluid that leaves it
        ///
        ///
        fn pushParticlesApart(&mut self, numIters: i32) {
            // // count particles per cell
            // this.numCellParticles.fill(0);
            self.numCellParticles.fill(0);

            // for (var i = 0; i < this.numParticles; i++) {
            // store all the particle positions and multiply by the grid spacing, make sure it's in the grid.
            for i in 0..self.numParticles {
                // var x = this.particlePos[2 * i];
                let x: f32 = self.particlePos[(2 * i) as usize];
                // var y = this.particlePos[2 * i + 1];
                let y: f32 = self.particlePos[(2 * i + 1) as usize];

                // var xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
                let xi = clamp(floorf(x) as i32, 1, self.pNumX - 2);
                // var yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
                let yi = clamp(floorf(y) as i32, 1, self.pNumY - 2);
                // var cellNr = xi * this.pNumY + yi;
                let celNr = xi * self.pNumY + yi;
                // this.numCellParticles[cellNr]++;
                self.numCellParticles[celNr as usize] += 1;
            }

            // // partial sums

            // var first = 0;
            let mut first = 0;

            // for (var i = 0; i < this.pNumCells; i++) {
            for i in 0..self.pNumCells {
                // first += this.numCellParticles[i];
                first += self.numCellParticles[i as usize];
                // this.firstCellParticle[i] = first;
                self.firstCellParticle[i as usize] = first;
            }
            // this.firstCellParticle[this.pNumCells] = first; // guard
            self.firstCellParticle[self.pNumCells as usize] = first;

            // // fill particles into cells

            // for (var i = 0; i < this.numParticles; i++) {
            for i in 0..self.numParticles {
                // var x = this.particlePos[2 * i];
                let x = self.particlePos[(2 * i) as usize];
                // var y = this.particlePos[2 * i + 1];
                let y = self.particlePos[(2 * i + 1) as usize];

                // var xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
                let xi = clamp(floorf(x * self.pInvSpacing) as i32, 1, self.pNumX - 2);
                // var yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
                let yi = clamp(floorf(y * self.pInvSpacing) as i32, 1, self.pNumY - 2);
                // var cellNr = xi * this.pNumY + yi;
                let cellNr = xi * self.pNumY + yi;
                // this.firstCellParticle[cellNr]--;
                self.firstCellParticle[cellNr as usize] -= 1;
                // this.cellParticleIds[this.firstCellParticle[cellNr]] = i;
                self.cellParticleIds[self.firstCellParticle[cellNr as usize] as usize] = i;
            }

            // // push particles apart

            // var minDist = 2.0 * this.particleRadius;
            let minDist = 2.0 * self.particleRadius;
            // var minDist2 = minDist * minDist;
            let minDist2 = minDist * minDist;

            // for (var iter = 0; iter < numIters; iter++) {
            for _i in 0..numIters {
                // for (var i = 0; i < this.numParticles; i++) {
                for i in 0..self.numParticles {
                    // var px = this.particlePos[2 * i];
                    let px = self.particlePos[(2 * i) as usize];
                    // var py = this.particlePos[2 * i + 1];
                    let py = self.particlePos[(2 * i + 1) as usize];

                    // var pxi = Math.floor(px * this.pInvSpacing);
                    let pxi = floorf(px * self.pInvSpacing) as i32;
                    // var pyi = Math.floor(py * this.pInvSpacing);
                    let pyi = floorf(py * self.pInvSpacing) as i32;
                    // var x0 = Math.max(pxi - 1, 0);
                    let x0 = (pxi - 1).max(0);
                    // var y0 = Math.max(pyi - 1, 0);
                    let y0 = (pyi - 1).max(0);
                    // var x1 = Math.min(pxi + 1, this.pNumX - 1);
                    let x1 = (pxi + 1).min(self.pNumX - 1);
                    // var y1 = Math.min(pyi + 1, this.pNumY - 1);
                    let y1 = (pyi + 1).min(self.pNumY - 1);

                    // for (var xi = x0; xi <= x1; xi++) {
                    for xi in x0..=x1 {
                        // for (var yi = y0; yi <= y1; yi++) {
                        for yi in y0..=y1 {
                            // var cellNr = xi * this.pNumY + yi;
                            let cellNr = xi * self.pNumY + yi;
                            // var first = this.firstCellParticle[cellNr];
                            let first = self.firstCellParticle[cellNr as usize];
                            // var last = this.firstCellParticle[cellNr + 1];
                            let last = self.firstCellParticle[(cellNr + 1) as usize];
                            // for (var j = first; j < last; j++) {
                            for j in first..last {
                                // var id = this.cellParticleIds[j];
                                let id = self.cellParticleIds[j as usize];
                                // if (id == i)
                                if id == i {
                                    continue;
                                    // continue;
                                }
                                if self.particlePos[(2 * id) as usize] == px
                                    && self.particlePos[(2 * id + 1) as usize] == py
                                {
                                    self.particlePos[2 * id as usize] += 0.001 * i as f32;
                                    self.particlePos[2 * id as usize + 1] += 0.001 * i as f32;
                                }
                                // var qx = this.particlePos[2 * id];
                                let qx = self.particlePos[(2 * id) as usize];
                                // var qy = this.particlePos[2 * id + 1];
                                let qy = self.particlePos[(2 * id + 1) as usize];

                                // var dx = qx - px;
                                let mut dx = qx - px;
                                // var dy = qy - py;
                                let mut dy = qy - py;
                                // var d2 = dx * dx + dy * dy;

                                let d2 = dx * dx + dy * dy;
                                // if (d2 > minDist2 || d2 == 0.0)
                                if d2 > minDist2 || d2 == 0f32 {
                                    // continue;
                                    continue;
                                }
                                // var d = Math.sqrt(d2);
                                let d = sqrtf(d2);
                                // var s = 0.5 * (minDist - d) / d;
                                let s = 0.5f32 * (minDist - d) / d;
                                // dx *= s;
                                dx *= s;
                                // dy *= s;
                                dy *= s;
                                // this.particlePos[2 * i] -= dx;
                                self.particlePos[(2 * i) as usize] -= dx;
                                // this.particlePos[2 * i + 1] -= dy;
                                self.particlePos[(2 * i + 1) as usize] -= dy;
                                // this.particlePos[2 * id] += dx;
                                self.particlePos[(2 * id) as usize] += dx;
                                // this.particlePos[2 * id + 1] += dy;
                                self.particlePos[(2 * id + 1) as usize] += dy;
                            }
                        }
                    }
                }
            }
        }

        // handleParticleCollisions(obstacleX, obstacleY, obstacleRadius) {
        fn handleParticleCollisions(
            &mut self,
            // obstacleX: f32,
            // obstacleY: f32,
            // obstacleRadius: f32,
            // scene: &Scene,
        ) {
            // var h = 1.0 / this.fInvSpacing;
            // var r = this.particleRadius;
            // var or = obstacleRadius;
            // let or = obstacleRadius;
            // var or2 = or * or;
            // let or2 = or * or;
            // var minDist = obstacleRadius + r;
            // let minDist = obstacleRadius + r;
            // var minDist2 = minDist * minDist;
            // let minDist2 = minDist * minDist;

            // var minX = h + r;
            let minX = 1.0;
            // var maxX = (this.fNumX - 1) * h - r;
            let maxX = 21.0;
            // var minY = h + r;
            let minY = 1.0;
            // var maxY = (this.fNumY - 1) * h - r;
            let maxY = 21.0;

            // for (var i = 0; i < this.numParticles; i++) {
            for i in 0..self.numParticles {
                // var x = this.particlePos[2 * i];
                let mut x = self.particlePos[(2 * i) as usize];
                // var y = this.particlePos[2 * i + 1];
                let mut y = self.particlePos[(2 * i + 1) as usize];

                // var dx = x - obstacleX;
                // let dx = x - obstacleX;
                // var dy = y - obstacleY;
                // let dy = y - obstacleY;
                // var d2 = dx * dx + dy * dy;
                // let d2 = dx + dy * dy;

                // // obstacle collision

                // if (d2 < minDist2) {
                // if d2 < minDist2 {
                // // var d = Math.sqrt(d2);
                // let d = sqrtf(d2);
                // // var s = (minDist - d) / d;
                // let s = (minDist - d) / d;
                // // x += dx * s;
                // x += dx * s;
                // // y += dy * s;
                // y += dy * s;

                // this.particleVel[2 * i] = scene.obstacleVelX;
                // self.particleVel[(2 * i) as usize] = scene.obstacleVelX;
                // this.particleVel[2 * i + 1] = scene.obstacleVelY;
                // self.particleVel[(2 * i + 1) as usize] = scene.obstacleVelY;
                // }

                // // wall collisions

                // if (x < minX) {
                if x < minX {
                    // x = minX;
                    x = minX;
                    // this.particleVel[2 * i] = 0.0;
                    self.particleVel[(2 * i) as usize] = 0.0;
                }
                // if (x > maxX) {
                if x > maxX {
                    // x = maxX;
                    x = maxX;
                    // this.particleVel[2 * i] = 0.0;
                    self.particleVel[(2 * i) as usize] = 0.0;
                }
                // if (y < minY) {
                if y < minY {
                    // y = minY;
                    y = minY;
                    // this.particleVel[2 * i + 1] = 0.0;
                    self.particleVel[(2 * i + 1) as usize] = 0.0;
                }
                // if (y > maxY) {
                if y > maxY {
                    // y = maxY;
                    y = maxY;
                    // this.particleVel[2 * i + 1] = 0.0;
                    self.particleVel[(2 * i + 1) as usize] = 0.0;
                }
                // this.particlePos[2 * i] = x;
                self.particlePos[(2 * i) as usize] = x;
                // this.particlePos[2 * i + 1] = y;
                self.particlePos[(2 * i + 1) as usize] = y;
            }
        }

        // updateParticleDensity() {
        fn updateParticleDensity(&mut self) {
            // var n = this.fNumY;
            let n = self.fNumY;
            // var h = this.h;
            let h = self.h;
            // var h1 = this.fInvSpacing;
            let h1 = self.fInvSpacing;
            // var h2 = 0.5 * h;
            let h2 = 0.5 * h;

            // var d = f.particleDensity;
            let mut d = self.particleDensity.clone();

            // d.fill(0.0);
            d.fill(0.0);

            // for (var i = 0; i < this.numParticles; i++) {
            for i in 0..self.numParticles {
                // var x = this.particlePos[2 * i];
                let mut x = self.particlePos[(2 * i) as usize];
                // var y = this.particlePos[2 * i + 1];
                let mut y = self.particlePos[(2 * i + 1) as usize];

                // x = clamp(x, h, (this.fNumX - 1) * h);
                x = clamp(x, h, (self.fNumX - 1.0) * h);
                // y = clamp(y, h, (this.fNumY - 1) * h);
                y = clamp(y, h, (self.fNumY - 1.0) * h);

                // var x0 = Math.floor((x - h2) * h1);
                let x0 = floorf((x - h2) * h1);
                // var tx = ((x - h2) - x0 * h) * h1;
                let tx = ((x - h2) - x0 * h) * h1;
                // var x1 = Math.min(x0 + 1, this.fNumX-2);
                let x1 = (x0 + 1.0).min(self.fNumX - 2.0);
                // var y0 = Math.floor((y-h2)*h1);
                let y0 = floorf((y - h2) * h1);
                // var ty = ((y - h2) - y0*h) * h1;
                let ty = ((y - h2) - y0 * h) * h1;
                // var y1 = Math.min(y0 + 1, this.fNumY-2);
                let y1 = (y0 + 1.0).min(self.fNumY - 2.0);

                // var sx = 1.0 - tx;
                let sx = 1.0 - tx;
                // var sy = 1.0 - ty;
                let sy = 1.0 - ty;

                // if (x0 < this.fNumX && y0 < this.fNumY) d[x0 * n + y0] += sx * sy;
                if x0 < self.fNumX && y0 < self.fNumY {
                    d[(x0 * n + y0) as usize] += sx * sy
                };
                // if (x1 < this.fNumX && y0 < this.fNumY) d[x1 * n + y0] += tx * sy;
                if x1 < self.fNumX && y0 < self.fNumY {
                    d[(x1 * n + y0) as usize] += tx * sy
                };
                // if (x1 < this.fNumX && y1 < this.fNumY) d[x1 * n + y1] += tx * ty;
                if x1 < self.fNumX && y1 < self.fNumY {
                    d[(x1 * n + y1) as usize] += tx * ty
                };
                // if (x0 < this.fNumX && y1 < this.fNumY) d[x0 * n + y1] += sx * ty;
                if x0 < self.fNumX && y1 < self.fNumY {
                    d[(x0 * n + y1) as usize] += sx * ty
                };
            }

            // if (this.particleRestDensity == 0.0) {
            if self.particleRestDensity == 0.0 {
                // var sum = 0.0;
                let mut sum = 0.0;
                // var numFluidCells = 0;
                let mut numFluidCells = 0;

                // for (var i = 0; i < this.fNumCells; i++) {
                for i in 0..self.fNumCells as usize {
                    // if (this.cellType[i] == FLUID_CELL) {
                    if self.cellType[i] == CellType::FLUID_CELL {
                        // sum += d[i];
                        sum += d[i];
                        // numFluidCells++;
                        numFluidCells += 1;
                    }
                }

                // if (numFluidCells > 0)
                if numFluidCells > 0 {
                    // this.particleRestDensity = sum / numFluidCells;
                    self.particleRestDensity = sum / numFluidCells as f32;
                }
            }

            // // for (var xi = 1; xi < this.fNumX; xi++) {

            // // for (var yi = 1; yi < this.fNumY; yi++) {
            // // var cellNr = xi * n + yi;
            // // if (this.cellType[cellNr] != FLUID_CELL)
            // // continue;
            // // var hx = this.h;
            // // var hy = this.h;

            // // if (this.cellType[(xi - 1) * n + yi] == SOLID_CELL || this.cellType[(xi + 1) * n + yi] == SOLID_CELL)
            // // hx -= this.particleRadius;
            // // if (this.cellType[xi * n + yi - 1] == SOLID_CELL || this.cellType[xi * n + yi + 1] == SOLID_CELL)
            // // hy -= this.particleRadius;

            // // var scale = this.h * this.h / (hx * hy)
            // // d[cellNr] *= scale;
            // // }
            // // }
        }

        // transferVelocities(toGrid, flipRatio){
        /// transfer the particle velocities to the grid or vice versa
        fn transferVelocities(&mut self, toGrid: bool, flipRatio: f32) {
            // var n = this.fNumY;
            // clone of the number of cells in the y direction (fNumY)
            let n = self.fNumY;
            // var h = this.h;
            // clone of the cell spacing (self.h)
            let h = self.h;
            // var h1 = this.fInvSpacing;
            // the inverse of the cell spacing (1/h)
            let h1 = self.fInvSpacing;
            // var h2 = 0.5 * h;
            // half the size of a cell
            let h2 = 0.5 * h;

            // clone cell positions and velocities into buffers and clear the values
            // for each cell, check if s = 0.0 and call it a solid cell if it is and air cell if not.
            // for each particle store the x and y position and
            // if (toGrid) {
            // set all the cell types, self.s = 0.0 => solid, if any particles in cell => fluid, else => air
            // (also stores previous positions and velocities in prevU/prevV and clears du/dv/u/v)
            if toGrid {
                // this.prevU.set(this.u);
                self.prevU = self.u.clone();
                // this.prevV.set(this.v);
                self.prevV = self.v.clone();

                // this.du.fill(0.0);
                self.du.fill(0.0);
                // this.dv.fill(0.0);
                self.dv.fill(0.0);
                // this.u.fill(0.0);
                self.u.fill(0.0);
                // this.v.fill(0.0);
                self.v.fill(0.0);

                // for (var i = 0; i < this.fNumCells; i++)
                // if s = 0.0 make it a solid cell, otherwise an air cell
                for i in 0..self.fNumCells as usize {
                    // this.cellType[i] = this.s[i] == 0.0 ? SOLID_CELL : AIR_CELL;
                    self.cellType[i] = if self.s[i] == 0.0 {
                        CellType::SOLID_CELL
                    } else {
                        CellType::AIR_CELL
                    };
                }

                // for (var i = 0; i < this.numParticles; i++) {
                // for each particle, get the cell that it's in and if that's an air cell make it a fluid cell
                for i in 0..self.numParticles {
                    // var x = this.particlePos[2 * i];
                    let x = self.particlePos[(2 * i) as usize];
                    // var y = this.particlePos[2 * i + 1];
                    let y = self.particlePos[(2 * i + 1) as usize];
                    // var xi = clamp(Math.floor(x * h1), 0, this.fNumX - 1);
                    let xi = clamp(floorf(x * h1), 1.0, self.fNumX - 2.0);
                    // var yi = clamp(Math.floor(y * h1), 0, this.fNumY - 1);
                    let yi = clamp(floorf(y * h1), 1.0, self.fNumY - 2.0);
                    // var cellNr = xi * n + yi;
                    let cellNr = xi * n + yi; // use the cell coordinates to get the ID of the cell it is in.
                    // if (this.cellType[cellNr] == AIR_CELL)
                    if self.cellType[cellNr as usize] == CellType::AIR_CELL {
                        // this.cellType[cellNr] = FLUID_CELL;
                        self.cellType[cellNr as usize] = CellType::FLUID_CELL;
                    }
                }
            }

            // for (var component = 0; component < 2; component++) {

            for component in 0..2 {
                // runs twice
                // var dx = component == 0 ? 0.0 : h2;
                let dx = if component == 0 { 0.0 } else { h2 }; //dx is 0 if component =0, and half the cell spacing if not
                // var dy = component == 0 ? h2 : 0.0;
                let dy = if component == 0 { h2 } else { 0.0 }; // dy is half the cell spacing if component=0 and 0 if not.

                // var f = component == 0 ? this.u : this.v;
                // on component=0, f is a clone of position (self.u), and on component=1 its a clone of velocity self.v
                let mut f = if component == 0 {
                    self.u.clone()
                } else {
                    self.v.clone()
                };
                // var prevF = component == 0 ? this.prevU : this.prevV;
                // on component = 0 prevF is a clone of prevU, and when component=1 it's a clone of prevV
                let prevF = if component == 0 {
                    self.prevU.clone()
                } else {
                    self.prevV.clone()
                };
                // var d = component == 0 ? this.du : this.dv;
                // on component = 0, d is a clone of du (the change in position per time) and for component=1 it's dv (change in velocity over time)
                let mut d = if component == 0 {
                    self.du.clone()
                } else {
                    self.dv.clone()
                };

                // for (var i = 0; i < this.numParticles; i++) {
                // so here there are 2 possibilities:
                // either component=0, dx=0  , dy=h/2, f=u.clone, prevF=prevU.clone, d=du.clone
                //     or component=1, dx=h/2, dy=0  , f=v.clone, prevF=prevV.clone, d=dv.clone
                //
                // for each particle
                for i in 0..self.numParticles {
                    // store the particle position in x and y
                    // var x = this.particlePos[2 * i];
                    let mut x = self.particlePos[(2 * i) as usize];
                    // var y = this.particlePos[2 * i + 1];
                    let mut y = self.particlePos[(2 * i + 1) as usize];

                    // clamp the positions between the cell spacing and the cell spacing times the grid size
                    // x = clamp(x, h, (this.fNumX - 1) * h);
                    x = clamp(x, h, (self.fNumX - 1.0) * h);
                    // y = clamp(y, h, (this.fNumY - 1) * h);
                    y = clamp(y, h, (self.fNumY - 1.0) * h);

                    //Xp = (x-dx)
                    //Xcell = floor(xp/h) -> x0 the x grid location
                    //DeltaX/h = (Xp-Xcell*h)/h -> tx
                    // x1 is x0 + 1 clamped to grid

                    // var x0 = Math.min(Math.floor((x - dx) * h1), this.fNumX - 2);
                    let x0 = floorf((x - dx) * h1).min(self.fNumX - 2.0);
                    // var tx = ((x - dx) - x0 * h) * h1;
                    let tx = ((x - dx) - x0 * h) * h1;
                    // var x1 = Math.min(x0 + 1, this.fNumX-2);
                    let x1 = (x0 + 1.0).min(self.fNumX - 2.0);

                    //Yp = (y-dy)
                    //Ycell = floor(yp/h) -> y0 the y grid location
                    //DeltaY/h = (Yp-Ycell*h)/h -> ty
                    // y1 is y0 + 1 clamped to grid

                    // var y0 = Math.min(Math.floor((y-dy)*h1), this.fNumY-2);
                    let y0 = floorf((y - dy) * h1).min(self.fNumY - 2.0);
                    // var ty = ((y - dy) - y0*h) * h1;
                    let ty = ((y - dy) - y0 * h) * h1;
                    // var y1 = Math.min(y0 + 1, this.fNumY-2);
                    let y1 = (y0 + 1.0).min(self.fNumY - 2.0);

                    //sx = 1-DeltaX/h
                    //sy = 1-DeltaY/h
                    // var sx = 1.0 - tx;
                    let sx = 1.0 - tx;
                    // var sy = 1.0 - ty;
                    let sy = 1.0 - ty;

                    // w1 = (1-DeltaX/h)(1-DeltaY/h) = sx * sy
                    // w2 = (DeltaX/h)(1-DeltaY/h) = tx * sy
                    // w3 = (DeltaX/h)(DeltaY/h) = tx * ty
                    // w4 = (1-DeltaX/h)(DeltaY/h) = sx * ty

                    // var d0 = sx*sy;
                    // w1
                    let d0 = sx * sy;
                    // var d1 = tx*sy;
                    // w2
                    let d1 = tx * sy;
                    // var d2 = tx*ty;
                    // w3
                    let d2 = tx * ty;
                    // var d3 = sx*ty;
                    // w4
                    let d3 = sx * ty;

                    // the x grid location times the number of cells in the y direction + the y grid location
                    // var nr0 = x0*n + y0;
                    let nr0 = x0 * n + y0; // q1, location of bottom left corner
                    // var nr1 = x1*n + y0;
                    let nr1 = x1 * n + y0; // q2, location of bottom right corner
                    // var nr2 = x1*n + y1;
                    let nr2 = x1 * n + y1; // q3, location of top right corner
                    // var nr3 = x0*n + y1;
                    let nr3 = x0 * n + y1; // q4, location of top left corner

                    // if (toGrid) {
                    if toGrid {
                        // var pv = this.particleVel[2 * i + component];
                        let pv = self.particleVel[(2 * i + component) as usize];
                        // f[nr0] += pv * d0; d[nr0] += d0;
                        f[nr0 as usize] += pv * d0;
                        d[nr0 as usize] += d0;
                        // f[nr1] += pv * d1; d[nr1] += d1;
                        f[nr1 as usize] += pv * d1;
                        d[nr1 as usize] += d1;
                        // f[nr2] += pv * d2; d[nr2] += d2;
                        f[nr2 as usize] += pv * d2;
                        d[nr2 as usize] += d2;
                        // f[nr3] += pv * d3; d[nr3] += d3;
                        f[nr3 as usize] += pv * d3;
                        d[nr3 as usize] += d3;
                    }
                    // else {
                    else {
                        // var offset = component == 0 ? n : 1;
                        let offset = if component == 0 { n } else { 1.0 };
                        // var valid0 = this.cellType[nr0] != AIR_CELL || this.cellType[nr0 - offset] != AIR_CELL ? 1.0 : 0.0;
                        let valid0 = if self.cellType[nr0 as usize] != CellType::AIR_CELL
                            || self.cellType[(nr0 - offset) as usize] != CellType::AIR_CELL
                        {
                            1.0
                        } else {
                            0.0
                        };
                        // var valid1 = this.cellType[nr1] != AIR_CELL || this.cellType[nr1 - offset] != AIR_CELL ? 1.0 : 0.0;
                        let valid1 = if self.cellType[nr1 as usize] != CellType::AIR_CELL
                            || self.cellType[(nr1 - offset) as usize] != CellType::AIR_CELL
                        {
                            1.0
                        } else {
                            0.0
                        };
                        // var valid2 = this.cellType[nr2] != AIR_CELL || this.cellType[nr2 - offset] != AIR_CELL ? 1.0 : 0.0;
                        let valid2 = if self.cellType[nr2 as usize] != CellType::AIR_CELL
                            || self.cellType[(nr2 - offset) as usize] != CellType::AIR_CELL
                        {
                            1.0
                        } else {
                            0.0
                        };
                        // var valid3 = this.cellType[nr3] != AIR_CELL || this.cellType[nr3 - offset] != AIR_CELL ? 1.0 : 0.0;
                        let valid3 = if self.cellType[nr3 as usize] != CellType::AIR_CELL
                            || self.cellType[(nr3 - offset) as usize] != CellType::AIR_CELL
                        {
                            1.0
                        } else {
                            0.0
                        };

                        // var v = this.particleVel[2 * i + component];
                        let v = self.particleVel[(2 * i + component) as usize];
                        // var d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;
                        let d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                        // if (d > 0.0) {
                        if d > 0.0 {
                            // var picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
                            let picV = (valid0 * d0 * f[nr0 as usize]
                                + valid1 * d1 * f[nr1 as usize]
                                + valid2 * d2 * f[nr2 as usize]
                                + valid3 * d3 * f[nr3 as usize])
                                / d;
                            // var corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1])
                            // + valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / d;
                            let corr = (valid0 * d0 * (f[nr0 as usize] - prevF[nr0 as usize])
                                + valid1 * d1 * (f[nr1 as usize] - prevF[nr1 as usize])
                                + valid2 * d2 * (f[nr2 as usize] - prevF[nr2 as usize])
                                + valid3 * d3 * (f[nr3 as usize] - prevF[nr3 as usize]))
                                / d;
                            // var flipV = v + corr;
                            let flipV = v + corr;

                            // this.particleVel[2 * i + component] = (1.0 - flipRatio) * picV + flipRatio * flipV;
                            self.particleVel[(2 * i + component) as usize] =
                                (1.0 - flipRatio) * picV + flipRatio * flipV;
                        }
                    }
                }

                // if (toGrid) {
                if toGrid {
                    // for (var i = 0; i < f.length; i++) {
                    for i in 0..f.len() {
                        // if (d[i] > 0.0)
                        if d[i] > 0.0 {
                            // f[i] /= d[i];
                            f[i] /= d[i];
                        }
                    }

                    // // restore solid cells

                    // for (var i = 0; i < this.fNumX; i++) {
                    for i in 0..self.fNumX as usize {
                        // for (var j = 0; j < this.fNumY; j++) {
                        for j in 0..self.fNumY as usize {
                            // var solid = this.cellType[i * n + j] == SOLID_CELL;
                            let solid = self.cellType[i * n as usize + j] == CellType::SOLID_CELL;
                            // if (solid || (i > 0 && this.cellType[(i - 1) * n + j] == SOLID_CELL))
                            if solid
                                || (i > 0
                                    && self.cellType[(i - 1) * n as usize + j]
                                        == CellType::SOLID_CELL)
                            {
                                // this.u[i * n + j] = this.prevU[i * n + j];
                                self.u[i * n as usize + j] = self.prevU[i * n as usize + j];
                            }
                            // if (solid || (j > 0 && this.cellType[i * n + j - 1] == SOLID_CELL))
                            if solid
                                || (j > 0
                                    && self.cellType[i * n as usize + j - 1]
                                        == CellType::SOLID_CELL)
                            {
                                // this.v[i * n + j] = this.prevV[i * n + j];
                                self.v[i * n as usize + j] = self.prevV[i * n as usize + j];
                            }
                        }
                    }
                }
            }
        }

        // solveIncompressibility(numIters, dt, overRelaxation, compensateDrift = true) {
        // make the grid velocities incompressible
        fn solveIncompressibility(
            &mut self,
            numIters: i32,
            dt: f32,
            overRelaxation: f32,
            compensateDrift: bool,
        ) {
            // this.p.fill(0.0);
            self.p.fill(0.0);
            // this.prevU.set(this.u);
            self.prevU = self.u.clone();
            // this.prevV.set(this.v);
            self.prevV = self.v.clone();

            // var n = this.fNumY;
            let n = self.fNumY;
            // var cp = this.density * this.h / dt;
            let cp = self.density * self.h / dt;

            // for (var i = 0; i < this.fNumCells; i++) {
            for i in 0..self.fNumCells as usize {
                // var u = this.u[i];
                let _u = self.u[i];
                // var v = this.v[i];
                let _v = self.v[i];
            }

            // for (var iter = 0; iter < numIters; iter++) {
            for _iter in 0..numIters {
                // for (var i = 1; i < this.fNumX-1; i++) {
                for i in 1..self.fNumX as usize {
                    // for (var j = 1; j < this.fNumY-1; j++) {
                    for j in 1..self.fNumY as usize {
                        // if (this.cellType[i*n + j] != FLUID_CELL)
                        if self.cellType[i * n as usize + j] != CellType::FLUID_CELL {
                            // continue;
                            continue;
                        }

                        // var center = i * n + j;
                        let center = i * n as usize + j;
                        // var left = (i - 1) * n + j;
                        let left = (i - 1) * n as usize + j;
                        // var right = (i + 1) * n + j;
                        let right = (i + 1) * n as usize + j;
                        // var bottom = i * n + j - 1;
                        let bottom = i * n as usize + j - 1;
                        // var top = i * n + j + 1;
                        let top = i * n as usize + j + 1;

                        // var s = this.s[center];
                        let _s = self.s[center];
                        // var sx0 = this.s[left];
                        let sx0 = self.s[left];
                        // var sx1 = this.s[right];
                        let sx1 = self.s[right];
                        // var sy0 = this.s[bottom];
                        let sy0 = self.s[bottom];
                        // var sy1 = this.s[top];
                        let sy1 = self.s[top];
                        // var s = sx0 + sx1 + sy0 + sy1;
                        let s = sx0 + sx1 + sy0 + sy1;
                        // if (s == 0.0)
                        if s == 0.0 {
                            // continue;
                            continue;
                        }

                        // var div = this.u[right] - this.u[center] +
                        // this.v[top] - this.v[center];
                        let mut div = self.u[right] - self.u[center] + self.v[top] - self.v[center];

                        // if (this.particleRestDensity > 0.0 && compensateDrift) {
                        if self.particleRestDensity > 0.0 && compensateDrift {
                            // var k = 1.0;
                            let k = 1.0;
                            // var compression = this.particleDensity[i*n + j] - this.particleRestDensity;
                            let compression =
                                self.particleDensity[i * n as usize + j] - self.particleRestDensity;
                            // if (compression > 0.0)
                            if compression > 0.0 {
                                // div = div - k * compression;
                                div = div - k * compression;
                            }
                        }

                        // var p = -div / s;
                        let mut p = -div / s;
                        // p *= overRelaxation;
                        p *= overRelaxation;
                        // this.p[center] += cp * p;
                        self.p[center] += cp * p;

                        // this.u[center] -= sx0 * p;
                        self.u[center] -= sx0 * p;
                        // this.u[right] += sx1 * p;
                        self.u[right] += sx1 * p;
                        // this.v[center] -= sy0 * p;
                        self.u[center] += sy0 * p;
                        // this.v[top] += sy1 * p;
                        self.v[top] += sy1 * p;
                    }
                }
            }
        }

        // simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters, overRelaxation, compensateDrift, separateParticles, obstacleX, abstacleY, obstacleRadius) {
        pub fn simulate(
            &mut self,
            dt: f32,
            xGravity: f32,
            yGravity: f32,
            flipRatio: f32,
            numPressureIters: i32,
            numParticleIters: i32,
            overRelaxation: f32,
            compensateDrift: bool,
            separateParticles: bool,
        ) {
            // var scene =

            // var numSubSteps = 1;
            let numSubSteps = 1.0;
            // var sdt = dt / numSubSteps;
            let sdt = dt / numSubSteps;

            // for (var step = 0; step < numSubSteps; step++) {
            for _i in 0..numSubSteps as usize {
                // this.integrateParticles(sdt, gravity);
                self.integrateParticles(sdt, yGravity, xGravity);
                // if (separateParticles)
                self.pushParticlesApart(numParticleIters);
                self.handleParticleCollisions();
                // // this.handleParticleCollisions(obstacleX, abstacleY, obstacleRadius)
                // // this.transferVelocities(true);
                self.transferVelocities(true, 1.9);
                // // this.updateParticleDensity();
                self.updateParticleDensity();
                // // this.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
                self.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
                // // this.transferVelocities(false, flipRatio);
                self.transferVelocities(false, flipRatio);
                self.showParticles();
            }

            // }
        }
    }

    pub struct Scene {
        xGravity: f32,
        yGravity: f32,
        dt: f32,
        flipRatio: f32,
        numPressureIters: i32,
        numParticleIters: i32,
        frameNr: i32,
        overRelaxation: f32,
        compensateDrift: bool,
        separateParticles: bool,
        _paused: bool,
        pub fluid: FlipFluid,
    }

    impl Scene {
        pub fn setupScene(particles: i32) -> Scene {
            // gravity : -9.81,
            let xGravity = 0.0;
            let yGravity = 0.0;
            // // gravity : 0.0,
            // // dt : 1.0 / 120.0,
            // let dt = 1.0 / 120.0;
            // flipRatio : 0.9,
            let flipRatio = 0.9;
            // // numPressureIters : 100,
            // let numPressureIters = 100;
            // // numParticleIters : 2,
            // let numParticleIters = 2;
            // frameNr : 0,
            let frameNr = 0;
            // overRelaxation : 1.9,
            let _overRelaxation = 1.9;
            // compensateDrift : true,
            let compensateDrift = true;
            // separateParticles : true,
            let separateParticles = true;
            // paused: true,
            let _paused = false;

            // scene.overRelaxation = 1.9;
            let overRelaxation = 1.9;
            // scene.dt = 1.0 / 60.0;
            let dt = 1.0 / 60.0;
            // scene.numPressureIters = 50;
            let numPressureIters = 5;
            // scene.numParticleIters = 2;
            let numParticleIters = 3;

            // var res = 100;
            let res = 23.0;
            // var tankHeight = 1.0 * simHeight;
            let tankHeight = 1.0 * simHeight;
            // var tankWidth = 1.0 * simWidth;
            let tankWidth = 1.0 * simWidth;
            // var h = tankHeight / res;
            let h = tankHeight / res;
            // var density = 1000.0;
            let density = 1000.0;

            // var relWaterHeight = 0.8
            let relWaterHeight = 0.8;
            // var relWaterWidth = 0.6
            let relWaterWidth = 0.6;

            // // dam break

            // // compute number of particles

            // var r = 0.3 * h; // particle radius w.r.t. cell size
            let r = 0.5 * h;
            // var dx = 2.0 * r;
            let dx = 2.0 * r;
            // var dy = Math.sqrt(3.0) / 2.0 * dx;
            let dy = sqrtf(3.0) / 2.0 * dx;

            // var numX = Math.floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
            let numX = floorf((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
            // var numY = Math.floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
            let numY = floorf((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
            // var maxParticles = numX * numY;
            let maxParticles = (numX * numY) as i32;

            // // create fluid

            // flash.blocking_write(reset_count_location, &[10u8]);
            // f = scene.fluid = new FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);
            let mut fluid = FlipFluid::new(density, tankWidth, tankHeight, h, r, maxParticles);

            fluid.numParticles = particles;

            Scene {
                xGravity,
                yGravity,
                dt,
                flipRatio,
                numPressureIters,
                numParticleIters,
                frameNr,
                overRelaxation,
                compensateDrift,
                separateParticles,
                _paused,
                fluid,
            }
            // // create particles

            // f.numParticles = numX * numY;
            // var p = 0;
            // for (var i = 0; i < numX; i++) {
            // for (var j = 0; j < numY; j++) {
            // f.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
            // f.particlePos[p++] = h + r + dy * j
            // }
            // }

            // // setup grid cells for tank

            // var n = f.fNumY;

            // for (var i = 0; i < f.fNumX; i++) {
            // for (var j = 0; j < f.fNumY; j++) {
            // var s = 1.0; // fluid
            // if (i == 0 || i == f.fNumX-1 || j == 0)
            // s = 0.0; // solid
            // f.s[i*n + j] = s
            // }
            // }

            // setObstacle(3.0, 2.0, true);
        }
        pub fn simulate(&mut self) {
            self.fluid.cellType.fill(CellType::AIR_CELL);
            self.fluid.simulate(
                self.dt,
                self.xGravity,
                self.yGravity,
                self.flipRatio,
                self.numPressureIters,
                self.numParticleIters,
                self.overRelaxation,
                self.compensateDrift,
                self.separateParticles,
            );
            self.frameNr += 1;
        }
        pub fn set_gravity(&mut self, accel_measurment: [f32; 2]) {
            self.xGravity = accel_measurment[0];
            self.yGravity = accel_measurment[1];
        }
        pub fn get_output(&mut self) -> [[bool; 21]; 21] {
            let mut output_frame: [[bool; 21]; 21] = [[false; 21]; 21];
            for i in 1..22 {
                for j in 1..22 {
                    if self.fluid.cellType[i * 23 + j] == CellType::FLUID_CELL {
                        output_frame[i - 1][j - 1] = true;
                    }
                }
            }
            output_frame
        }
    }

    // function clamp(x, min, max) {
    fn clamp<T: PartialOrd>(x: T, min: T, max: T) -> T {
        // if (x < min)
        if x < min {
            // return min;
            min
        }
        // else if (x > max)
        else if x > max {
            // return max;
            max
        }
        // else
        else {
            // return x;
            x
        }
    }
}
