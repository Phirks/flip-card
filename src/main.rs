#![no_std]
#![no_main]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused)]
#![allow(non_upper_case_globals)]
use core::num;
use embassy_executor::Spawner;
use embassy_futures::select::select;
use embassy_rp::bind_interrupts;
use embassy_rp::gpio;
use embassy_rp::i2c;
use embassy_rp::i2c::Config as I2C_Config;
use embassy_rp::i2c::InterruptHandler as I2C_InterruptHandler;
use embassy_rp::i2c::{AbortReason, Error};
use embassy_rp::pac;
use embassy_rp::pac::common::W;
use embassy_rp::pac::rosc::regs::Count;
use embassy_rp::peripherals::I2C1;
use embassy_rp::peripherals::PIO0;
use embassy_rp::pio::program::pio_asm;
use embassy_rp::pio::{Config, InterruptHandler, Pio, ShiftConfig, ShiftDirection};
use embassy_rp::watchdog::*;
use embassy_sync::watch;
use embassy_time::Timer;
use embedded_hal_async::i2c::I2c;
use fixed::traits::ToFixed;
use fixed_macro::types::U56F8;
use gpio::{Level, Output};
use heapless::Vec;
use libm::{exp, floor, floorf, sin, sqrtf};
// use {defmt_rtt as _, panic_probe as _};

static max_particles_setting: usize = 400;
static number_of_vertical_cells_setting: usize = 23;
static number_of_horizontal_cells_setting: usize = 23;

static max_particles_x2_setting: usize = max_particles_setting * 3;
static max_particles_x3_setting: usize = max_particles_setting * 3;
static number_of_cells_setting: usize =
    number_of_vertical_cells_setting * number_of_horizontal_cells_setting;
static number_of_cells_x2_setting: usize = number_of_cells_setting * 2;
static number_of_cells_x3_setting: usize = number_of_cells_setting * 3;
static canvas_height: f32 = number_of_vertical_cells_setting as f32;
static canvas_width: f32 = number_of_horizontal_cells_setting as f32;
// var canvas = document.getElementById("myCanvas");
// var gl = canvas.getContext("webgl");
// canvas.width = window.innerWidth - 20;
// canvas.height = window.innerHeight - 20;

// canvas.focus();

// var simHeight = 3.0;
static simHeight: f32 = 3.0;
// var cScale = canvas.height / simHeight;
static cScale: f32 = canvas_height / simHeight;
// var simWidth = canvas.width / cScale;
static simWidth: f32 = canvas_width / cScale;

// var U_FIELD = 0;
static U_FIELD: i32 = 0;
// var V_FIELD = 1;
static V_FIELD: i32 = 1;

static BIT21: u32 = 0b00000000001000000000000000000000; // pinout 19 1
static BIT20: u32 = 0b00000000000100000000000000000000; // pinout 19 2
static BIT19: u32 = 0b00000000000010000000000000000000; // pinout 19 3
static BIT18: u32 = 0b00000000000001000000000000000000; // pinout 18 4
static BIT17: u32 = 0b00000000000000100000000000000000; // pinout 17 5
static BIT16: u32 = 0b00000000000000010000000000000000; // pinout 16 6
static BIT15: u32 = 0b00000000000000001000000000000000; // pinout 15 7
static BIT14: u32 = 0b00000000000000000100000000000000; // pinout 14 8
static BIT13: u32 = 0b00000000000000000010000000000000; // pinout 13 9
static BIT12: u32 = 0b00000000000000000001000000000000; // pinout 12 10
static BIT11: u32 = 0b00000000000000000000100000000000; // pinout 11 11
static BIT10: u32 = 0b00000000000000000000010000000000; // pinout 10 12
static BIT9: u32 = 0b00000000000000000000001000000000; // pinout 9 13
static BIT8: u32 = 0b00000000000000000000000100000000; // pinout 8 14
static BIT7: u32 = 0b00000000000000000000000010000000; // pinout 7 15
static BIT6: u32 = 0b00000000000000000000000001000000; // pinout 6 16
static BIT5: u32 = 0b00000000000000000000000000100000; // pinout 5 17
static BIT4: u32 = 0b00000000000000000000000000010000; // pinout 4 18
static BIT3: u32 = 0b00000000000000000000000000001000; // pinout 3 19
static BIT2: u32 = 0b00000000000000000000000000000100; // pinout 2 20
static BIT1: u32 = 0b00000000000000000000000000000010; // pinout 1 21
static BIT0: u32 = 0b00000000000000000000000000000001; // pinout 1 21

static BITS: [u32; 22] = [
    BIT21, BIT20, BIT19, BIT18, BIT17, BIT16, BIT15, BIT14, BIT13, BIT12, BIT11, BIT10, BIT9, BIT8,
    BIT7, BIT6, BIT5, BIT4, BIT3, BIT2, BIT1, BIT0,
];

static LOOKUP_TABLE: [[usize; 21]; 21] = [
    [
        416, 481, 373, 435, 327, 389, 281, 343, 235, 297, 189, 251, 143, 205, 097, 159, 051, 113,
        005, 067, 025,
    ],
    [
        459, 372, 479, 329, 433, 283, 387, 237, 341, 191, 295, 145, 249, 099, 203, 053, 157, 007,
        111, 027, 068,
    ],
    [
        370, 457, 328, 477, 285, 431, 239, 385, 193, 339, 147, 293, 101, 247, 055, 201, 009, 155,
        029, 112, 071,
    ],
    [
        413, 326, 455, 284, 475, 241, 429, 195, 383, 149, 337, 103, 291, 057, 245, 011, 199, 031,
        156, 073, 114,
    ],
    [
        324, 411, 282, 453, 240, 473, 197, 427, 151, 381, 105, 335, 059, 289, 013, 243, 033, 200,
        075, 158, 117,
    ],
    [
        367, 280, 409, 238, 451, 196, 471, 153, 425, 107, 379, 061, 333, 015, 287, 035, 244, 077,
        202, 119, 160,
    ],
    [
        278, 365, 236, 407, 194, 449, 152, 469, 109, 423, 063, 377, 017, 331, 037, 288, 079, 246,
        121, 204, 163,
    ],
    [
        321, 234, 363, 192, 405, 150, 447, 108, 467, 065, 421, 019, 375, 039, 332, 081, 290, 123,
        248, 165, 206,
    ],
    [
        232, 319, 190, 361, 148, 403, 106, 445, 064, 465, 021, 419, 041, 376, 083, 334, 125, 292,
        167, 250, 209,
    ],
    [
        275, 188, 317, 146, 359, 104, 401, 062, 443, 020, 463, 043, 420, 085, 378, 127, 336, 169,
        294, 211, 252,
    ],
    [
        186, 273, 144, 315, 102, 357, 060, 399, 018, 441, 042, 464, 087, 422, 129, 380, 171, 338,
        213, 296, 255,
    ],
    [
        229, 142, 271, 100, 313, 058, 355, 016, 397, 040, 442, 086, 466, 131, 424, 173, 382, 215,
        340, 257, 298,
    ],
    [
        140, 227, 098, 269, 056, 311, 014, 353, 038, 398, 084, 444, 130, 468, 175, 426, 217, 384,
        259, 342, 301,
    ],
    [
        183, 096, 225, 054, 267, 012, 309, 036, 354, 082, 400, 128, 446, 174, 470, 219, 428, 261,
        386, 303, 344,
    ],
    [
        094, 181, 052, 223, 010, 265, 034, 310, 080, 356, 126, 402, 172, 448, 218, 472, 263, 430,
        305, 388, 347,
    ],
    [
        137, 050, 179, 008, 221, 032, 266, 078, 312, 124, 358, 170, 404, 216, 450, 262, 474, 307,
        432, 349, 390,
    ],
    [
        048, 135, 006, 177, 030, 222, 076, 268, 122, 314, 168, 360, 214, 406, 260, 452, 306, 476,
        351, 434, 393,
    ],
    [
        091, 004, 133, 028, 178, 074, 224, 120, 270, 166, 316, 212, 362, 258, 408, 304, 454, 350,
        478, 395, 436,
    ],
    [
        002, 089, 026, 134, 072, 180, 118, 226, 164, 272, 210, 318, 256, 364, 302, 410, 348, 456,
        394, 480, 439,
    ],
    [
        045, 024, 090, 070, 136, 116, 182, 162, 228, 208, 274, 254, 320, 300, 366, 346, 412, 392,
        458, 438, 482,
    ],
    [
        023, 046, 069, 092, 115, 138, 161, 184, 207, 230, 253, 276, 299, 322, 345, 368, 391, 414,
        437, 460, 461,
    ],
];

// static QR_CODE: [[u8; 21]; 21] = [
//     [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
//     [0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
//     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
// ];

static QR_CODE: [[u8; 21]; 21] = [
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
    [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ],
];

#[derive(PartialEq)]
enum CellType {
    AIR_CELL,
    FLUID_CELL,
    SOLID_CELL,
}

struct FlipFluid {
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
    u: Vec<f32, number_of_cells_setting>,

    /// an array of the cell velocities 2*i is the vertical velocity, 2*i+1 is the horizontal velocity
    v: Vec<f32, number_of_cells_setting>,

    /// this is an array of the amount the position will change in one time step 2*i is the
    /// y coordinate change and 2*i+1 is the x coordinate change
    du: Vec<f32, number_of_cells_setting>,

    /// this is an array of the amount the velocity will change in one time step 2*i is the
    /// vertical velocity change and 2*i+1 is the horizontal velocity change.
    dv: Vec<f32, number_of_cells_setting>,

    /// this is an array of the previous iteration's particle position 2*i are y coordinates
    /// and 2*i+1 are x coordinates
    prevU: Vec<f32, number_of_cells_setting>,

    /// this is an array of the previous iteration's particle velocity.  2*i are y vertical
    /// velocities and 2*i+1 are horizontal velocities.
    prevV: Vec<f32, number_of_cells_setting>,

    p: Vec<f32, number_of_cells_setting>,
    s: Vec<f32, number_of_cells_setting>,
    cellType: Vec<CellType, number_of_cells_setting>,
    cellColor: Vec<f32, number_of_cells_x3_setting>,

    /// the max number of particles, used to size the arrays
    maxParticles: i32,

    /// the particle position, in meters, indexes 2*i are vertical, 2*i+1 are horizontal
    particlePos: Vec<f32, max_particles_x2_setting>,

    /// the particle color, in groups of 3 for RGB
    particleColor: Vec<f32, max_particles_x3_setting>,

    /// the particle velocity, in m/s, indexes 2*i are vertical, 2*i+1 are horizontal
    particleVel: Vec<f32, max_particles_x2_setting>,

    particleDensity: Vec<f32, number_of_cells_setting>,
    particleRestDensity: f32,
    particleRadius: f32,

    pInvSpacing: f32,
    pNumX: i32,
    pNumY: i32,
    pNumCells: i32,

    numCellParticles: Vec<i32, 30>,
    firstCellParticle: Vec<i32, 32>,
    cellParticleIds: Vec<i32, max_particles_setting>,

    /// the total number of particles in the system
    numParticles: i32,
}

impl FlipFluid {
    fn new(
        density: f32,
        width: f32,
        height: f32,
        spacing: f32,
        particleRadius: f32,
        maxParticles: i32,
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
        let fNumCells = (fNumX * fNumY);

        // this.u = new Float32Array(this.fNumCells);
        let u = Vec::new();
        // this.v = new Float32Array(this.fNumCells);
        let v = Vec::new();
        // this.du = new Float32Array(this.fNumCells);
        let du = Vec::new();
        // this.dv = new Float32Array(this.fNumCells);
        let dv = Vec::new();
        // this.prevU = new Float32Array(this.fNumCells);
        let prevU = Vec::new();
        // this.prevV = new Float32Array(this.fNumCells);
        let prevV = Vec::new();
        // this.p = new Float32Array(this.fNumCells);
        let p = Vec::new();
        // this.s = new Float32Array(this.fNumCells);
        let s = Vec::new();
        // this.cellType = new Int32Array(this.fNumCells);
        let cellType = Vec::new();
        // this.cellColor = new Float32Array(3 * this.fNumCells);
        let cellColor = Vec::new();
        // this.particleRadius = particleRadius;
        // this.maxParticles = maxParticles;
        //maxParticles
        // this.particlePos = new Float32Array(2 * this.maxParticles);
        let particlePos = Vec::new();
        // this.particleColor = new Float32Array(3 * this.maxParticles);
        let mut particleColor = Vec::new();
        // for (var i = 0; i < this.maxParticles; i++)
        // this.particleColor[3 * i + 2] = 1.0;
        for i in 0..maxParticles {
            particleColor[(3 * i + 2) as usize] = 1.0;
        }
        // this.particleVel = new Float32Array(2 * this.maxParticles);
        let particleVel = Vec::new();
        // this.particleDensity = new Float32Array(this.fNumCells);
        let particleDensity = Vec::new();
        // this.particleRestDensity = 0.0;
        let particleRestDensity = 0.0f32;

        //particleRadius
        // this.pInvSpacing = 1.0 / (2.2 * particleRadius);
        let pInvSpacing = 1.0 / (2.2 * particleRadius);
        // this.pNumX = Math.floor(width * this.pInvSpacing) + 1;
        let pNumX = floorf(width * pInvSpacing) as i32;
        // this.pNumY = Math.floor(height * this.pInvSpacing) + 1;
        let pNumY = floorf(height * pInvSpacing) as i32;
        // this.pNumCells = this.pNumX * this.pNumY;
        let pNumCells = pNumX * pNumY;
        // this.numCellParticles = new Int32Array(this.pNumCells);
        let numCellParticles = Vec::new();
        // this.firstCellParticle = new Int32Array(this.pNumCells + 1);
        let cellParticleIds = Vec::new();
        // this.cellParticleIds = new Int32Array(maxParticles);
        let firstCellParticle = Vec::new();
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
            cellColor,
            maxParticles,
            particlePos,
            particleColor,
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
    fn integrateParticles(&mut self, dt: f32, gravity: f32) {
        // for (var i = 0; i < this.numParticles; i++) {
        for i in 0..self.numParticles {
            // this.particleVel[2 * i + 1] += dt * gravity;
            self.particleVel[(2 * i + 1) as usize] += dt * gravity;
            // this.particlePos[2 * i] += this.particleVel[2 * i] * dt;
            self.particlePos[(2 * i) as usize] += self.particleVel[(2 * i) as usize] * dt;
            // this.particlePos[2 * i + 1] += this.particleVel[2 * i + 1] * dt;
            self.particlePos[(2 * i + 1) as usize] += self.particleVel[(2 * i + 1) as usize] * dt;
        }
    }
    // pushParticlesApart(numIters)
    /// store the particle positions in x and y and
    /// make incompressible by making sure the amount of fluid that enters a cell is equal to the fluid that leaves it
    ///
    ///
    fn pushParticlesApart(&mut self, numIters: i32) {
        // var colorDiffusionCoeff = 0.001;
        let colorDiffusionCoeff: f32 = 0.001;
        // // count particles per cell
        // this.numCellParticles.fill(0);
        self.numCellParticles.fill(0);

        // for (var i = 0; i < this.numParticles; i++) {
        for i in 0..self.numParticles {
            // var x = this.particlePos[2 * i];
            let x: f32 = self.particlePos[(2 * i) as usize];
            // var y = this.particlePos[2 * i + 1];
            let y: f32 = self.particlePos[(2 * i + 1) as usize];

            // var xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
            let xi = clamp(floorf(x * self.pInvSpacing) as i32, 0, self.pNumX - 1);
            // var yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
            let yi = clamp(floorf(y * self.pInvSpacing) as i32, 0, self.pNumY - 1);
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
            let xi = clamp(floorf(x * self.pInvSpacing) as i32, 0, self.pNumX);
            // var yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
            let yi = clamp(floorf(y * self.pInvSpacing) as i32, 0, self.pNumY);
            // var cellNr = xi * this.pNumY + yi;
            let cellNr = xi * self.pNumY + yi;
            // this.firstCellParticle[cellNr]--;
            self.firstCellParticle[cellNr as usize];
            // this.cellParticleIds[this.firstCellParticle[cellNr]] = i;
            self.cellParticleIds[self.firstCellParticle[cellNr as usize] as usize] = i;
        }

        // // push particles apart

        // var minDist = 2.0 * this.particleRadius;
        let minDist = 2.0 * self.particleRadius;
        // var minDist2 = minDist * minDist;
        let minDist2 = minDist * minDist;

        // for (var iter = 0; iter < numIters; iter++) {
        for i in 0..numIters {
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
                for xi in x0..x1 {
                    // for (var yi = y0; yi <= y1; yi++) {
                    for yi in y0..y1 {
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

                            // // diffuse colors

                            // for (var k = 0; k < 3; k++) {
                            for k in 0..3 {
                                // var color0 = this.particleColor[3 * i + k];
                                let color0 = self.particleColor[(3 * i + k) as usize];
                                // var color1 = this.particleColor[3 * id + k];
                                let color1 = self.particleColor[(3 * id + k) as usize];
                                // var color = (color0 + color1) * 0.5;
                                let color = (color0 + color1) * 0.5f32;
                                // this.particleColor[3 * i + k] = color0 + (color - color0) * colorDiffusionCoeff;
                                self.particleColor[(3 * i + k) as usize] =
                                    color0 + (color - color0) * colorDiffusionCoeff;
                                // this.particleColor[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
                                self.particleColor[(3 * id + k) as usize] =
                                    color1 + (color - color1) * colorDiffusionCoeff;
                            }
                        }
                    }
                }
            }
        }
    }

    // handleParticleCollisions(obstacleX, obstacleY, obstacleRadius) {
    fn handleParticleCollisions(
        &mut self,
        obstacleX: f32,
        obstacleY: f32,
        obstacleRadius: f32,
        scene: &Scene,
    ) {
        // var h = 1.0 / this.fInvSpacing;
        let h = 1.0 / self.fInvSpacing;
        // var r = this.particleRadius;
        let r = self.particleRadius;
        // var or = obstacleRadius;
        let or = obstacleRadius;
        // var or2 = or * or;
        let or2 = or * or;
        // var minDist = obstacleRadius + r;
        let minDist = obstacleRadius + r;
        // var minDist2 = minDist * minDist;
        let minDist2 = minDist * minDist;

        // var minX = h + r;
        let minX = h + r;
        // var maxX = (this.fNumX - 1) * h - r;
        let maxX = (self.fNumX - 1.0) * h - r;
        // var minY = h + r;
        let minY = h + r;
        // var maxY = (this.fNumY - 1) * h - r;
        let maxY = (self.fNumY - 1.0) * h - r;

        // for (var i = 0; i < this.numParticles; i++) {
        for i in 0..self.numParticles {
            // var x = this.particlePos[2 * i];
            let mut x = self.particlePos[(2 * i) as usize];
            // var y = this.particlePos[2 * i + 1];
            let mut y = self.particlePos[(2 * i + 1) as usize];

            // var dx = x - obstacleX;
            let dx = x - obstacleX;
            // var dy = y - obstacleY;
            let dy = y - obstacleY;
            // var d2 = dx * dx + dy * dy;
            let d2 = dx + dy * dy;

            // // obstacle collision

            // if (d2 < minDist2) {
            if d2 < minDist2 {
                // // var d = Math.sqrt(d2);
                let d = sqrtf(d2);
                // // var s = (minDist - d) / d;
                let s = (minDist - d) / d;
                // // x += dx * s;
                x += dx * s;
                // // y += dy * s;
                y += dy * s;

                // this.particleVel[2 * i] = scene.obstacleVelX;
                self.particleVel[(2 * i) as usize] = scene.obstacleVelX;
                // this.particleVel[2 * i + 1] = scene.obstacleVelY;
                self.particleVel[(2 * i + 1) as usize] = scene.obstacleVelY;
            }

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
        let n = self.fNumY;
        // var h = this.h;
        let h = self.h;
        // var h1 = this.fInvSpacing;
        let h1 = self.fInvSpacing;
        // var h2 = 0.5 * h;
        let h2 = 0.5 * h;

        // clone cell positions and velocities into buffers and clear the values
        // for each cell, check if s = 0.0 and call it a solid cell if it is and air cell if not.
        // for each particle store the x and y position and
        // if (toGrid) {
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
            for i in 0..self.fNumCells as usize {
                // this.cellType[i] = this.s[i] == 0.0 ? SOLID_CELL : AIR_CELL;
                self.cellType[i] = if self.s[i] == 0.0 {
                    CellType::SOLID_CELL
                } else {
                    CellType::AIR_CELL
                };
            }

            // for (var i = 0; i < this.numParticles; i++) {
            for i in 0..self.numParticles {
                // var x = this.particlePos[2 * i];
                let x = self.particlePos[(2 * i) as usize];
                // var y = this.particlePos[2 * i + 1];
                let y = self.particlePos[(2 * i + 1) as usize];
                // var xi = clamp(Math.floor(x * h1), 0, this.fNumX - 1);
                let xi = clamp(floorf(x * h1), 0.0, self.fNumX - 1.0);
                // var yi = clamp(Math.floor(y * h1), 0, this.fNumY - 1);
                let yi = clamp(floorf(y * h1), 0.0, self.fNumY - 1.0);
                // var cellNr = xi * n + yi;
                let cellNr = xi * n + yi;
                // if (this.cellType[cellNr] == AIR_CELL)
                if self.cellType[cellNr as usize] == CellType::AIR_CELL {
                    // this.cellType[cellNr] = FLUID_CELL;
                    self.cellType[cellNr as usize] = CellType::FLUID_CELL;
                }
            }
        }

        // for (var component = 0; component < 2; component++) {
        for component in 0..2 {
            // var dx = component == 0 ? 0.0 : h2;
            let dx = if component == 0 { 0.0 } else { h2 };
            // var dy = component == 0 ? h2 : 0.0;
            let dy = if component == 0 { h2 } else { 0.0 };

            // var f = component == 0 ? this.u : this.v;
            let mut f = if component == 0 {
                self.u.clone()
            } else {
                self.v.clone()
            };
            // var prevF = component == 0 ? this.prevU : this.prevV;
            let prevF = if component == 0 {
                self.prevU.clone()
            } else {
                self.prevV.clone()
            };
            // var d = component == 0 ? this.du : this.dv;
            let mut d = if component == 0 {
                self.du.clone()
            } else {
                self.dv.clone()
            };

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

                // var x0 = Math.min(Math.floor((x - dx) * h1), this.fNumX - 2);
                let x0 = ((x - dx) * h1).min(self.fNumX - 2.0);
                // var tx = ((x - dx) - x0 * h) * h1;
                let tx = ((x - dx) - x0 * h) * h1;
                // var x1 = Math.min(x0 + 1, this.fNumX-2);
                let x1 = (x0 + 1.0).min(self.fNumX - 2.0);
                // var y0 = Math.min(Math.floor((y-dy)*h1), this.fNumY-2);
                let y0 = floorf((y - dy) * h1).min(self.fNumY - 2.0);
                // var ty = ((y - dy) - y0*h) * h1;
                let ty = ((y - dy) - y0 * h) * h1;
                // var y1 = Math.min(y0 + 1, this.fNumY-2);
                let y1 = (y0 + 1.0).min(self.fNumY - 2.0);

                // var sx = 1.0 - tx;
                let sx = 1.0 - tx;
                // var sy = 1.0 - ty;
                let sy = 1.0 - ty;

                // var d0 = sx*sy;
                let d0 = sx * sy;
                // var d1 = tx*sy;
                let d1 = tx * sy;
                // var d2 = tx*ty;
                let d2 = tx * ty;
                // var d3 = sx*ty;
                let d3 = sx * ty;

                // var nr0 = x0*n + y0;
                let nr0 = x0 * n + y0;
                // var nr1 = x1*n + y0;
                let nr1 = x1 * n + y0;
                // var nr2 = x1*n + y1;
                let nr2 = x1 * n + y1;
                // var nr3 = x0*n + y1;
                let nr3 = x0 * n + y1;

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
                                && self.cellType[((i - 1) * n as usize + j)]
                                    == CellType::SOLID_CELL)
                        {
                            // this.u[i * n + j] = this.prevU[i * n + j];
                            self.u[i * n as usize + j] = self.prevU[i * n as usize + j];
                        }
                        // if (solid || (j > 0 && this.cellType[i * n + j - 1] == SOLID_CELL))
                        if solid
                            || (j > 0
                                && self.cellType[i * n as usize + j - 1] == CellType::SOLID_CELL)
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
            let u = self.u[i];
            // var v = this.v[i];
            let v = self.v[i];
        }

        // for (var iter = 0; iter < numIters; iter++) {
        for iter in 0..numIters {
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
                    let s = self.s[center];
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

    // updateParticleColors() {
    fn updateParticleColors(&mut self) {
        // // for (var i = 0; i < this.numParticles; i++) {
        // // this.particleColor[3 * i] *= 0.99;
        // // this.particleColor[3 * i + 1] *= 0.99
        // // this.particleColor[3 * i + 2] =
        // // clamp(this.particleColor[3 * i + 2] + 0.001, 0.0, 1.0)
        // // }

        // // return;

        // var h1 = this.fInvSpacing;
        let h1 = self.fInvSpacing;

        // for (var i = 0; i < this.numParticles; i++) {
        for i in 0..self.numParticles as usize {
            // var s = 0.01;
            let s: f32 = 0.01;
            // this.particleColor[3 * i] = clamp(this.particleColor[3 * i] - s, 0.0, 1.0);
            self.particleColor[(3 * i)] = clamp(self.particleColor[(3 * i)] - s, 0.0, 1.0);
            // this.particleColor[3 * i + 1] = clamp(this.particleColor[3 * i + 1] - s, 0.0, 1.0);
            self.particleColor[(3 * i + 1)] = clamp(self.particleColor[3 * i + 1] - s, 0.0, 1.0);
            // this.particleColor[3 * i + 2] = clamp(this.particleColor[3 * i + 2] + s, 0.0, 1.0);
            self.particleColor[(3 * i + 1)] = clamp(self.particleColor[3 * i + 2] + s, 0.0, 1.0);

            // var x = this.particlePos[2 * i];
            let x = self.particlePos[2 * i];
            // var y = this.particlePos[2 * i + 1];
            let y = self.particlePos[2 * i + 1];
            // var xi = clamp(Math.floor(x * h1), 1, this.fNumX - 1);
            let xi = clamp(floorf(x * h1), 1.0, self.fNumX - 1.0);
            // var yi = clamp(Math.floor(y * h1), 1, this.fNumY - 1);
            let yi = clamp(floorf(y * h1), 1.0, self.fNumY - 1.0);
            // var cellNr = xi * this.fNumY + yi;
            let cellNr = xi * self.fNumY + yi;

            // var d0 = this.particleRestDensity;
            let d0 = self.particleRestDensity;

            // if (d0 > 0.0) {
            if d0 > 0.0 {
                // var relDensity = this.particleDensity[cellNr] / d0;
                let relDensity = self.particleDensity[cellNr as usize] / d0;
                // if (relDensity < 0.7) {
                if relDensity < 0.7 {
                    // var s = 0.8;
                    let s = 0.8;
                    // this.particleColor[3 * i] = s;
                    self.particleColor[3 * i] = s;
                    // this.particleColor[3 * i + 1] = s;
                    self.particleColor[3 * i + 1] = s;
                    // this.particleColor[3 * i + 2] = 1.0;
                    self.particleColor[3 * i + 2] = 1.0;
                }
            }
        }
    }

    // setSciColor(cellNr, val, minVal, maxVal) {
    fn setSciColor(&mut self, cellNr: usize, val: f32, minVal: f32, maxVal: f32) {
        // val = Math.min(Math.max(val, minVal), maxVal- 0.0001);
        let mut val = val.max(minVal).min(maxVal - 0.0001);
        // var d = maxVal - minVal;
        let d = maxVal - minVal;
        // val = d == 0.0 ? 0.5 : (val - minVal) / d;
        let val = if d == 0.0 { 0.5 } else { (val - minVal) / d };
        // var m = 0.25;
        let m: f32 = 0.25;
        // var num = Math.floor(val / m);
        let num = floorf(val / m);
        // var s = (val - num * m) / m;
        let s = (val - num * m) / m;
        // var r, g, b;
        let mut r: f32 = 0.0;
        let mut g: f32 = 0.0;
        let mut b: f32 = 0.0;

        // switch (num) {
        match num {
            // case 0 : r = 0.0; g = s; b = 1.0; break;
            0.0 => {
                r = 0.0;
                g = s;
                b = 1.0;
            }
            // case 1 : r = 0.0; g = 1.0; b = 1.0-s; break;
            1.0 => {
                r = 0.0;
                g = 1.0;
                b = 1.0 - s;
            }
            // case 2 : r = s; g = 1.0; b = 0.0; break;
            2.0 => {
                r = s;
                g = 1.0;
                b = 0.0;
            }
            // case 3 : r = 1.0; g = 1.0 - s; b = 0.0; break;
            3.0 => {
                r = 1.0;
                g = 1.0 - s;
                b = 0.0;
            }
            _ => {}
        }

        // this.cellColor[3 * cellNr] = r;
        self.cellColor[3 * cellNr as usize] = r;
        // this.cellColor[3 * cellNr + 1] = g;
        self.cellColor[3 * cellNr as usize + 1] = g;
        // this.cellColor[3 * cellNr + 2] = b;
        self.cellColor[3 * cellNr as usize + 2] = b;
    }

    // updateCellColors() {
    fn updateCellColors(&mut self) {
        // this.cellColor.fill(0.0);
        self.cellColor.fill(0.0);

        // for (var i = 0; i < this.fNumCells; i++) {
        for i in 0..self.fNumCells as usize {
            // if (this.cellType[i] == SOLID_CELL) {
            if self.cellType[i] == CellType::SOLID_CELL {
                // this.cellColor[3*i] = 0.5;
                self.cellColor[3 * i] = 0.5;
                // this.cellColor[3*i + 1] = 0.5;
                self.cellColor[3 * i + 1] = 0.5;
                // this.cellColor[3*i + 2] = 0.5;
                self.cellColor[3 * i + 2] = 0.5;
            }
            // else if (this.cellType[i] == FLUID_CELL) {
            else if self.cellType[i] == CellType::FLUID_CELL {
                // var d = this.particleDensity[i];
                let mut d = self.particleDensity[i];
                // if (this.particleRestDensity > 0.0)
                if self.particleRestDensity > 0.0 {
                    // d /= this.particleRestDensity;
                    d /= self.particleRestDensity;
                }
                // this.setSciColor(i, d, 0.0, 2.0);
                self.setSciColor(i, d, 0.0, 2.0);
            }
        }
    }

    // simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters, overRelaxation, compensateDrift, separateParticles, obstacleX, abstacleY, obstacleRadius) {
    fn simulate(
        &mut self,
        dt: f32,
        gravity: f32,
        flipRatio: f32,
        numPressureIters: i32,
        numParticleIters: i32,
        overRelaxation: f32,
        compensateDrift: bool,
        separateParticles: bool,
        obstacleX: f32,
        obstacleY: f32,
        obstacleRadius: f32,
    ) {
        // var scene =

        // var numSubSteps = 1;
        let numSubSteps = 1.0;
        // var sdt = dt / numSubSteps;
        let sdt = dt / numSubSteps;

        // for (var step = 0; step < numSubSteps; step++) {
        for i in 0..numSubSteps as usize {
            // this.integrateParticles(sdt, gravity);
            self.integrateParticles(sdt, gravity);
            // if (separateParticles)
            if separateParticles {
                // this.pushParticlesApart(numParticleIters);
                self.pushParticlesApart(numParticleIters);
            }
            // this.handleParticleCollisions(obstacleX, abstacleY, obstacleRadius)
            //           self.handleParticleCollisions(obstacleX, obstacleY, obstacleRadius, &scene);
            // this.transferVelocities(true);
            self.transferVelocities(true, 1.9);
            // this.updateParticleDensity();
            self.updateParticleDensity();
            // this.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
            self.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
            // this.transferVelocities(false, flipRatio);
            self.transferVelocities(false, flipRatio);
        }

        // this.updateParticleColors();
        self.updateParticleColors();
        // this.updateCellColors();
        self.updateCellColors();

        // }
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

enum FluidType {
    Water,
}

struct Scene {
    gravity: f32,
    dt: f32,
    flipRatio: f32,
    numPressureIters: i32,
    numParticleIters: i32,
    frameNr: i32,
    overRelaxation: f32,
    compensateDrift: bool,
    separateParticles: bool,
    obstacleX: f32,
    obstacleY: f32,
    obstacleRadius: f32,
    paused: bool,
    showObstacle: bool,
    obstacleVelX: f32,
    obstacleVelY: f32,
    showParticles: bool,
    showGrid: bool,
    fluid: FlipFluid,
}

impl Scene {
    fn setupScene() -> Scene {
        // gravity : -9.81,
        let gravity = -9.81;
        // // gravity : 0.0,
        // dt : 1.0 / 120.0,
        let dt = 1.0 / 120.0;
        // flipRatio : 0.9,
        let flipRatio = 0.9;
        // numPressureIters : 100,
        let numPressureIters = 100;
        // numParticleIters : 2,
        let numParticleIters = 2;
        // frameNr : 0,
        let frameNr = 0;
        // overRelaxation : 1.9,
        let overRelaxation = 1.9;
        // compensateDrift : true,
        let compensateDrift = true;
        // separateParticles : true,
        let separateParticles = true;
        // obstacleX : 0.0,
        let obstacleX = 0.0;
        // obstacleY : 0.0,
        let obstacleY = 0.0;
        // obstacleRadius: 0.15,
        let obstacleRadius = 0.15;
        // paused: true,
        let paused = true;
        // showObstacle: true,
        let showObstacle = true;
        // obstacleVelX: 0.0,
        let obstacleVelX = 0.0;
        // obstacleVelY: 0.0,
        let obstacleVelY = 0.0;
        // showParticles: true,
        let showParticles = true;
        // showGrid: false,
        let showGrid = false;
        // fluid: null

        // function setupScene() {
        // scene.obstacleRadius = 0.15;
        let obstacleRadius = 0.15;
        // scene.overRelaxation = 1.9;
        let overRelaxation = 1.9;
        // scene.dt = 1.0 / 60.0;
        let dt = 1.0 / 60.0;
        // scene.numPressureIters = 50;
        let numPressureIters = 50;
        // scene.numParticleIters = 2;
        let numParticleIters = 2;

        // var res = 100;
        let res = 100.0;
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
        let r = 0.3 * h;
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

        // f = scene.fluid = new FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);
        let fluid = FlipFluid::new(density, tankWidth, tankHeight, h, r, maxParticles);

        Scene {
            gravity,
            dt,
            flipRatio,
            numPressureIters,
            numParticleIters,
            frameNr,
            overRelaxation,
            compensateDrift,
            separateParticles,
            obstacleX,
            obstacleY,
            obstacleRadius,
            paused,
            showObstacle,
            obstacleVelX,
            obstacleVelY,
            showParticles,
            showGrid,
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
}

struct Screen {
    yx_grid: [[bool; 21]; 21],
    index_grid: [bool; 484],
    out_array: [u32; 484],
}
impl Screen {
    fn fill_index(&mut self) {
        self.index_grid = [false; 484];
        for y in 0..21 {
            for x in 0..21 {
                if self.yx_grid[y][x] {
                    self.index_grid[LOOKUP_TABLE[y][x]] = true
                }
            }
        }
    }
    fn one_pixel(&mut self, (y, x): (usize, usize)) {
        self.index_grid = [false; 484];
        self.index_grid[LOOKUP_TABLE[y][x]] = true;
    }
    fn make_output(&mut self) {
        self.out_array = [0; 484];
        let mut index = 0;

        for i in 0..BITS.len() {
            self.out_array[index] = BITS[i];
            index += 1;
            for j in 0..BITS.len() {
                if i != j {
                    if self.index_grid[index] {
                        self.out_array[index] = BITS[i] | BITS[j];
                    }
                    index += 1;
                }
            }
        }
    }
    fn test_item(&mut self, num: usize) {
        //self.out_array = [0; 484];
        let mut index = 0;
        for i in 0..BITS.len() {
            self.out_array[index] = BITS[i];
            index += 1;
            for j in 0..BITS.len() {
                if i != j {
                    if index == num {
                        self.out_array[index] = BITS[i] | BITS[j];
                    }
                    index += 1;
                }
            }
        }
    }
    fn disp_num(&mut self, num: usize, yloc: usize, xloc: usize) {
        let ymin = yloc;
        let ymax = yloc + 7;
        let xmin = xloc;
        let xmax = xloc + 5;

        for y in ymin..ymax {
            for x in xmin..xmax {
                self.yx_grid[y][x] = false;
            }
        }
        match num {
            0 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                //self.yx_grid[ymin +3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            1 => {
                // self.yx_grid[ymin + 1][xmin + 1] = true;
                // self.yx_grid[ymin + 2][xmin + 1] = true;
                // self.yx_grid[ymin + 3][xmin + 1] = true;
                // self.yx_grid[ymin + 4][xmin + 1] = true;
                // self.yx_grid[ymin + 5][xmin + 1] = true;

                // self.yx_grid[ymin + 1][xmin + 2] = true;
                // self.yx_grid[ymin +3][xmin + 2] = true;
                // self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            2 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                // self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                // self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            3 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                // self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                // self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            4 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                // self.yx_grid[ymin + 4][xmin + 1] = true;
                // self.yx_grid[ymin + 5][xmin + 1] = true;

                // self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                // self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            5 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                // self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                // self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            6 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                // self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                // self.yx_grid[ymin + 1][xmin + 3] = true;
                // self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            7 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                // self.yx_grid[ymin + 2][xmin + 1] = true;
                // self.yx_grid[ymin + 3][xmin + 1] = true;
                // self.yx_grid[ymin + 4][xmin + 1] = true;
                // self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                // self.yx_grid[ymin +3][xmin + 2] = true;
                // self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            8 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            9 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                // self.yx_grid[ymin + 4][xmin + 1] = true;
                // self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                // self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            10 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                // self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            11 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                //self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                //self.yx_grid[ymin + 1][xmin + 3] = true;
                //self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            12 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                //self.yx_grid[ymin +3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                //self.yx_grid[ymin + 2][xmin + 3] = true;
                //self.yx_grid[ymin + 3][xmin + 3] = true;
                //self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            13 => {
                //self.yx_grid[ymin + 1][xmin + 1] = true;
                //self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                //self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            14 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                // self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                // self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            15 => {
                self.yx_grid[ymin + 1][xmin + 1] = true;
                self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                // self.yx_grid[ymin + 5][xmin + 2] = true;

                self.yx_grid[ymin + 1][xmin + 3] = true;
                // self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                // self.yx_grid[ymin + 4][xmin + 3] = true;
                // self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            16 => {
                // self.yx_grid[ymin + 1][xmin + 1] = true;
                // self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                // self.yx_grid[ymin + 4][xmin + 1] = true;
                self.yx_grid[ymin + 5][xmin + 1] = true;

                // self.yx_grid[ymin + 1][xmin + 2] = true;
                // self.yx_grid[ymin +3][xmin + 2] = true;
                // self.yx_grid[ymin + 5][xmin + 2] = true;
                self.yx_grid[ymin + 4][xmin + 2] = true;

                // self.yx_grid[ymin + 1][xmin + 3] = true;
                // self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                // self.yx_grid[ymin + 4][xmin + 3] = true;
                self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            17 => {
                // self.yx_grid[ymin + 1][xmin + 1] = true;
                // self.yx_grid[ymin + 2][xmin + 1] = true;
                self.yx_grid[ymin + 3][xmin + 1] = true;
                // self.yx_grid[ymin + 4][xmin + 1] = true;
                // self.yx_grid[ymin + 5][xmin + 1] = true;

                // self.yx_grid[ymin + 1][xmin + 2] = true;
                self.yx_grid[ymin + 3][xmin + 2] = true;
                // self.yx_grid[ymin + 5][xmin + 2] = true;

                // self.yx_grid[ymin + 1][xmin + 3] = true;
                // self.yx_grid[ymin + 2][xmin + 3] = true;
                self.yx_grid[ymin + 3][xmin + 3] = true;
                // self.yx_grid[ymin + 4][xmin + 3] = true;
                // self.yx_grid[ymin + 5][xmin + 3] = true;
            }
            _ => {}
        }
    }
}

// Program metadata for `picotool info`.
// This isn't needed, but it's recomended to have these minimal entries.
#[unsafe(link_section = ".bi_entries")]
#[used]
pub static PICOTOOL_ENTRIES: [embassy_rp::binary_info::EntryAddr; 4] = [
    embassy_rp::binary_info::rp_program_name!(c"Blinky Example"),
    embassy_rp::binary_info::rp_program_description!(
        c"This example tests the RP Pico on board LED, connected to gpio 25"
    ),
    embassy_rp::binary_info::rp_cargo_version!(),
    embassy_rp::binary_info::rp_program_build_attribute!(),
];

bind_interrupts!(struct Irqs {  //sets up the IRQ for the PIO, probably to pull data out of the RX FIFO
    PIO0_IRQ_0 => InterruptHandler<PIO0>;
    I2C1_IRQ => I2C_InterruptHandler<I2C1>;
});

#[embassy_executor::main]
async fn main(_spawner: Spawner) {
    let p = embassy_rp::init(Default::default());
    // let mut watchdog = Watchdog::new(p.WATCHDOG);
    let mut enable_accel = Output::new(p.PIN_29, Level::Low);
    enable_accel.set_high();
    let pio = p.PIO0; //this PIO object is one of the 4 PIO peripherals
    let Pio {
        //why is this capitalized?
        mut common,
        sm0: mut sm,
        ..
    } = Pio::new(pio, Irqs); //contains the pio peripheral

    let i2c_config = embassy_rp::i2c::Config::default();
    let mut i2c = embassy_rp::i2c::I2c::new_async(p.I2C1, p.PIN_23, p.PIN_22, Irqs, i2c_config);

    let addr: u8 = 0x18;
    // set ctrl register
    i2c.write(addr, &[0x20, 0x57]).await;

    let prg = pio_asm!(
        //this program  is for moving data from TX to RX
        "out pins, 32", //change the state of pins, 0 is off, 1 is on.
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
        "out pindirs, 32",
    );

    let mut cfg = Config::default(); // how does this know it's a PIO config?

    let p0 = common.make_pio_pin(p.PIN_0);
    let p1 = common.make_pio_pin(p.PIN_1);
    let p2 = common.make_pio_pin(p.PIN_2);
    let p3 = common.make_pio_pin(p.PIN_3);
    let p4 = common.make_pio_pin(p.PIN_4);
    let p5 = common.make_pio_pin(p.PIN_5);
    let p6 = common.make_pio_pin(p.PIN_6);
    let p7 = common.make_pio_pin(p.PIN_7);
    let p8 = common.make_pio_pin(p.PIN_8);
    let p9 = common.make_pio_pin(p.PIN_9);
    let p10 = common.make_pio_pin(p.PIN_10);
    let p11 = common.make_pio_pin(p.PIN_11);
    let p12 = common.make_pio_pin(p.PIN_12);
    let p13 = common.make_pio_pin(p.PIN_13);
    let p14 = common.make_pio_pin(p.PIN_14);
    let p15 = common.make_pio_pin(p.PIN_15);
    let p16 = common.make_pio_pin(p.PIN_16);
    let p17 = common.make_pio_pin(p.PIN_17);
    let p18 = common.make_pio_pin(p.PIN_18);
    let p19 = common.make_pio_pin(p.PIN_19);
    let p20 = common.make_pio_pin(p.PIN_20);
    let p21 = common.make_pio_pin(p.PIN_21);

    cfg.set_out_pins(&[
        &p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8, &p9, &p10, &p11, &p12, &p13, &p14, &p15, &p16,
        &p17, &p18, &p19, &p20, &p21,
    ]);
    cfg.use_program(&common.load_program(&prg.program), &[]); // load the program instructions
    cfg.clock_divider = (U56F8!(125_000_000) / 300_000).to_fixed(); //set the speed, I thing its 12,500Hz
    cfg.shift_in = ShiftConfig {
        //config for shift in data from the left into the TX
        auto_fill: true, //autofills the tx buffer when 32 bits shifted out
        threshold: 32,
        direction: ShiftDirection::Right,
    };
    cfg.shift_out = ShiftConfig {
        //shift data out from the right into the RX
        auto_fill: true,
        threshold: 32,
        direction: ShiftDirection::Left,
    };
    sm.set_config(&cfg);
    sm.set_enable(true);

    let mut dma_out_ref = p.DMA_CH0;

    let mut screen = Screen {
        yx_grid: [[false; 21]; 21],
        index_grid: [false; 484],
        out_array: [0; 484],
    };
    //screen.yx_grid[2][20] = false;
    for y in 0..21 {
        for x in 0..21 {
            if QR_CODE[y][x] == 1 {
                screen.yx_grid[y][x] = true;
            } else {
                screen.yx_grid[y][x] = false;
            }
        }
    }

    let tx = sm.tx(); // this is where the state machine interface is
    let mut xh: [u8; 1] = [0];
    let mut yh: [u8; 1] = [0];
    let mut zh: [u8; 1] = [0];
    let mut x_val: i8 = 0;
    let mut y_val: i8 = 0;
    let mut z_val: i8 = 0;
    let mut x_negative: usize = 18;
    let mut y_negative: usize = 18;
    let mut z_negative: usize = 18;
    let mut dot = Dot::new();
    loop {
        i2c.write_read(addr, &[0x29], &mut xh).await; // read accel data
        i2c.write_read(addr, &[0x2B], &mut yh).await; // read accel data
        // i2c.write_read(addr, &[0x2D], &mut zh).await; // read accel data

        x_val = xh[0] as i8;
        y_val = yh[0] as i8;
        // z_val = zh[0] as i8;

        // if x_val < 0 {
        //     x_val = -x_val;
        //     x_negative = 17;
        // } else {
        //     x_negative = 18;
        // }
        // if y_val < 0 {
        //     y_val = -y_val;
        //     y_negative = 17;
        // } else {
        //     y_negative = 18;
        // }
        // if z_val < 0 {
        //     z_val = -z_val;
        //     z_negative = 17;
        // } else {
        //     z_negative = 18;
        // }
        if y_val > 60 {
            embassy_rp::rom_data::reboot(0x0002, 1, 0x00, 0x00); // reboot to BOOTSEL
        }
        dot.accelerate((-y_val as f32) / 1000.0, (x_val as f32) / 1000.0);
        screen.one_pixel(dot.grid_location);
        // screen.fill_index();
        // screen.disp_num(x_negative as usize, 0, 0);
        // screen.disp_num((x_val / 100).unsigned_abs() as usize, 0, 4);
        // screen.disp_num((x_val / 10 % 10).unsigned_abs() as usize, 0, 8);
        // screen.disp_num((x_val % 10).unsigned_abs() as usize, 0, 12);
        // screen.disp_num(y_negative as usize, 7, 0);
        // screen.disp_num((y_val / 100).unsigned_abs() as usize, 7, 4);
        // screen.disp_num((y_val / 10 % 10).unsigned_abs() as usize, 7, 8);
        // screen.disp_num((y_val % 10).unsigned_abs() as usize, 7, 12);
        // screen.disp_num(z_negative as usize, 14, 0);
        // screen.disp_num((z_val / 100).unsigned_abs() as usize, 14, 4);
        // screen.disp_num((z_val / 10 % 10).unsigned_abs() as usize, 14, 8);
        // screen.disp_num((z_val % 10).unsigned_abs() as usize, 14, 12);

        screen.make_output();
        for _ in 0..3 {
            tx.dma_push(dma_out_ref.reborrow(), &screen.out_array, false)
                .await;
        }
    }
}

struct Dot {
    x: f32,
    y: f32,
    x_velocity: f32,
    y_velocity: f32,
    grid_location: (usize, usize),
}
impl Dot {
    fn new() -> Dot {
        Dot {
            x: 0.0,
            y: 0.0,
            y_velocity: 0.0,
            x_velocity: 0.0,
            grid_location: (0, 0),
        }
    }
    fn accelerate(&mut self, y_accel: f32, x_accel: f32) {
        self.y_velocity += y_accel;
        self.x_velocity += x_accel;
        self.x += self.x_velocity;
        self.y += self.y_velocity;
        self.y = clamp(self.y, 0.0, 20.0);
        self.x = clamp(self.x, 0.0, 20.0);
        self.grid_location = (floorf(self.y) as usize, floorf(self.x) as usize);
        if self.x == 0.0 || self.x == 20.0 {
            self.x_velocity = 0.0
        };
        if self.y == 0.0 || self.y == 20.0 {
            self.y_velocity = 0.0
        };
    }
}

#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! {
    loop {}
}
