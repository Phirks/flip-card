#![no_std]
#![no_main]

use core::default;

use embassy_executor::{Executor, Spawner};
use embassy_futures::select::{Either::First, select};
use embassy_rp::bind_interrupts;
use embassy_rp::gpio;
use embassy_rp::gpio::DormantWakeConfig;
use embassy_rp::gpio::Input;
use embassy_rp::gpio::Pull;
use heapless::Vec;
// use embassy_rp::gpio::{Input, Pull};
use embassy_rp::clocks::{ClockConfig, CoreVoltage, clk_sys_freq, core_voltage};
use embassy_rp::i2c;
use embassy_rp::i2c::I2c;
use embassy_rp::i2c::InterruptHandler as I2C_InterruptHandler;
use embassy_rp::multicore::{Stack, spawn_core1};
use embassy_rp::peripherals::*;
use embassy_rp::pio::program::pio_asm;
use embassy_rp::pio::{Config, InterruptHandler, Pio, ShiftConfig, ShiftDirection};
use embassy_rp::watchdog::*;
use embassy_sync::blocking_mutex::raw::CriticalSectionRawMutex;
use embassy_time::Duration;
use embassy_time::Timer;
use fixed::traits::ToFixed;
use fixed_macro::types::U56F8;
use fluid_sim::FluidSimulation::Scene;
use gpio::{Level, Output};

use static_cell::StaticCell;
// use {defmt_rtt as _, panic_probe as _};

static mut CORE1_STACK: Stack<262144> = Stack::new();
static EXECUTOR0: StaticCell<Executor> = StaticCell::new();
static EXECUTOR1: StaticCell<Executor> = StaticCell::new();
// static CHANNEL: Channel<CriticalSectionRawMutex, LedState, 1> = Channel::new();

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

static FRAME_DATA_SIGNAL: embassy_sync::signal::Signal<CriticalSectionRawMutex, [[bool; 21]; 21]> =
    embassy_sync::signal::Signal::new();

static ACCEL_DATA_SIGNAL: embassy_sync::signal::Signal<CriticalSectionRawMutex, [f32; 2]> =
    embassy_sync::signal::Signal::new();

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

struct Screen {
    yx_grid: [[bool; 21]; 21],
    yx_lock_grid: [[bool; 21]; 21],
    index_grid: [bool; 484],
    out_array: [u32; 506],
    out_array_2: Vec<u32, 1452>,
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
    fn post_process(&mut self) {
        let mut temp_yx_lock_grid = [[false; 21]; 21];
        // First pass: populate temp_yx_lock_grid
        for i in 0..21 {
            for j in 0..21 {
                if self.yx_grid[i][j] {
                    temp_yx_lock_grid[i][j] = 
                        (i == 20 || self.yx_grid[i + 1][j]) &&
                        (i == 0  || self.yx_grid[i - 1][j]) &&
                        (j == 20 || self.yx_grid[i][j + 1]) &&
                        (j == 0  || self.yx_grid[i][j - 1]);
                }
            }
        }
        
        // Second pass: update yx_grid based on yx_lock_grid
        for i in 0..21 {
            for j in 0..21 {
                if self.yx_lock_grid[i][j] {
                    self.yx_grid[i][j] = 
                        (i == 20 || self.yx_grid[i + 1][j]) &&
                        (i == 0  || self.yx_grid[i - 1][j]) &&
                        (j == 20 || self.yx_grid[i][j + 1]) &&
                        (j == 0  || self.yx_grid[i][j - 1]);
                }
            }
        }
        self.yx_lock_grid = temp_yx_lock_grid;
    }
    fn make_output(&mut self) {
        self.out_array = [0; 506];
        let mut index = 0;
        let mut counter = 0;

        for i in 0..BITS.len() {
            self.out_array[counter] = BITS[i];
            index += 1;
            counter += 1;
            for j in 0..BITS.len() {
                if i != j {
                    if self.index_grid[index] {
                        self.out_array[counter] = BITS[i] | BITS[j];
                    }
                    index += 1;
                    counter += 1;
                }
            }
            // counter += 1;
        }
    }
    fn make_output_2(&mut self) {
        self.out_array_2 = Vec::new();
        let mut index = 0;
        for i in 0..BITS.len() {
            index += 1;
            for j in 0..BITS.len() {
                if i != j {
                    if self.index_grid[index] {
                        self.out_array_2.push(BITS[i]).unwrap();
                        self.out_array_2.push(BITS[i] | BITS[j]).unwrap();
                        self.out_array_2.push(0).unwrap();
                    }
                    index += 1;
                }
            }
        }
        if self.index_grid[461] {
            self.out_array_2.push(BITS[21]).unwrap();
            self.out_array_2.push(BITS[21] | BITS[20]).unwrap();
            self.out_array_2.push(0).unwrap();
        }
    }
    fn _disp_num(&mut self, num: usize, yloc: usize, xloc: usize) {
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
    embassy_rp::binary_info::rp_program_name!(c"flip-card"),
    embassy_rp::binary_info::rp_program_description!(
        c"This is a business card that runs a fluid simulation on an LED grid"
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
    // define all the peripherals and pins
    let mut config =
        embassy_rp::config::Config::new(ClockConfig::system_freq(200_000_000).unwrap());
    config.clocks.core_voltage = CoreVoltage::V1_15;
    let p = embassy_rp::init(config);
    let mut i2c: I2c<'static, I2C1, i2c::Async> =
        embassy_rp::i2c::I2c::new_async(p.I2C1, p.PIN_23, p.PIN_22, Irqs, Default::default());
    let watchdog = Watchdog::new(p.WATCHDOG);
    let mut enable_accel = Output::new(p.PIN_29, Level::Low);
    let pio: embassy_rp::Peri<'static, PIO0> = p.PIO0; //this PIO object is one of the 4 PIO peripherals
    let dma_out_ref: embassy_rp::Peri<'static, DMA_CH0> = p.DMA_CH0;
    let mut accel_interrupt_1 = Input::new(p.PIN_24, Pull::Up);
    let mut accel_interrupt_2 = Input::new(p.PIN_25, Pull::Up);
    // Timer::after_millis(10).await;

    // start up the peripherals
    let sm = setup_pio_2(
        pio, p.PIN_0, p.PIN_1, p.PIN_2, p.PIN_3, p.PIN_4, p.PIN_5, p.PIN_6, p.PIN_7, p.PIN_8,
        p.PIN_9, p.PIN_10, p.PIN_11, p.PIN_12, p.PIN_13, p.PIN_14, p.PIN_15, p.PIN_16, p.PIN_17,
        p.PIN_18, p.PIN_19, p.PIN_20, p.PIN_21,
    );
    enable_accel.set_high();
    let screen = Screen {
        yx_grid: [[false; 21]; 21],
        yx_lock_grid: [[false; 21]; 21],
        index_grid: [false; 484],
        out_array: [0; 506],
        out_array_2: Vec::new(),
    };

    // spawn background tasks

    spawn_core1(
        p.CORE1,
        unsafe { &mut *core::ptr::addr_of_mut!(CORE1_STACK) },
        move || {
            let executor1 = EXECUTOR1.init(Executor::new());
            executor1.run(|spawner| spawner.spawn(simulation_update()).unwrap())
        },
    );

    let executor0 = EXECUTOR0.init(Executor::new());
    executor0.run(|spawner| {
        spawner
            .spawn(monitor_accelerometer(
                i2c,
                accel_interrupt_1,
                accel_interrupt_2,
                enable_accel,
            ))
            .unwrap();
        spawner
            .spawn(drive_screen(watchdog, screen, sm, dma_out_ref))
            .unwrap()
    });
}

#[embassy_executor::task(pool_size = 1)]
async fn drive_screen(
    mut watchdog: Watchdog,
    mut screen: Screen,
    mut sm: embassy_rp::pio::StateMachine<'static, PIO0, 0>,
    mut dma_out_ref: embassy_rp::Peri<'static, DMA_CH0>,
) {
    let tx = sm.tx();
    let mut frame_count: u128 = 0;
    FRAME_DATA_SIGNAL.wait().await;
    watchdog.start(Duration::from_millis(500));
    loop {
        Timer::after_millis(1).await;
        watchdog.feed();
        if FRAME_DATA_SIGNAL.signaled() {
            screen.yx_grid = FRAME_DATA_SIGNAL.wait().await;
            frame_count = 0;
        } else {
            frame_count += 1;
        }
        screen.post_process();
        screen.fill_index();
        // screen.make_output();
        screen.make_output_2();
        tx.dma_push(dma_out_ref.reborrow(), &screen.out_array_2, false)
            .await;
        if frame_count >= 30 {
            // flash.blocking_write(reset_count_location, &[10u8]);
            // embassy_rp::rom_data::reboot(0x0002, 1, 0x00, 0x01); // reboot to BOOTSEL
            Timer::after_millis(3000).await;
        }
    }
}

#[embassy_executor::task(pool_size = 1)]
async fn monitor_accelerometer(
    mut i2c: I2c<'static, I2C1, i2c::Async>,
    mut int1: Input<'static>,
    mut int2: Input<'static>,
    mut enable_accel: Output<'static>,
) {
    //setup accelerometer
    let addr: u8 = 0x18;
    let dt: f32 = 1.0 / 100.0; // 1/Hz,
    let accel_scale_max: f32 = 9.81 * 8.0; // m/s²
    let mut xh: [u8; 1] = [0];
    i2c.blocking_write(addr, &[0x24, 0x08]).unwrap(); // CTRL_REG5, reboot chip
    i2c.blocking_write(addr, &[0x24, 0x80]).unwrap(); // CTRL_REG1, turn on, enable xyz
    Timer::after_millis(100).await;
    i2c.blocking_write(addr, &[0x20, 0x77]).unwrap(); // CTRL_REG1, turn on, enable xyz
    i2c.blocking_write(addr, &[0x21, 0x00]).unwrap(); // CTRL_REG2, high pass filter off
    i2c.blocking_write(addr, &[0x22, 0x80]).unwrap(); // CTRL_REG3, click to INT1 pin
    i2c.blocking_write(addr, &[0x23, 0x20]).unwrap(); // CTRL_REG4, FS = 8G
    // i2c.blocking_write(addr, &[0x24, 0x08]).unwrap(); // CTRL_REG5, latch INT pin 1
    i2c.blocking_write(addr, &[0x32, 0x0A]).unwrap(); // INT1_THS, Threshold = 250 mg
    i2c.blocking_write(addr, &[0x33, 0x7F]).unwrap(); // INT1_DURATION, Duration = 0
    i2c.blocking_write(addr, &[0x38, 0x22]).unwrap(); // CLICK_CFG
    i2c.blocking_write(addr, &[0x39, 0x27]).unwrap(); // CLICK_SRC
    i2c.blocking_write(addr, &[0x3A, 0x40]).unwrap(); // CLICK_THS
    i2c.blocking_write(addr, &[0x3B, 0x20]).unwrap(); // TIME_LIMIT
    i2c.blocking_write(addr, &[0x3C, 0x2F]).unwrap(); // TIME_LATENCY
    i2c.blocking_write(addr, &[0x3D, 0x30]).unwrap(); // TIME WINDOW
    // i2c.blocking_write(addr, &[0x30, 0x02]).unwrap(); // INT1_CFG, enable xh and yh interrupts
    // int1.wait_for_high().await;
    enable_accel.set_low();
    let dormant_config = DormantWakeConfig {
        edge_high: false,
        edge_low: false,
        level_high: true,
        level_low: false,
    };
    let wake = int1.dormant_wake(dormant_config.clone());
    let wake2 = int2.dormant_wake(dormant_config.clone());
    embassy_rp::clocks::dormant_sleep();
    drop(wake);
    drop(wake2);

    enable_accel.set_high();
    // i2c.blocking_write(addr, &[0x24, 0x80]).unwrap(); // CTRL_REG5, reboot chip
    i2c.blocking_write(addr, &[0x20, 0x57]).unwrap(); // CTRL_REG1, turn on, enable xy
    i2c.blocking_write(addr, &[0x23, 0x20]).unwrap(); // CTRL_REG1, turn on, enable xy
    //read accelerometer
    let mut yh: [u8; 1] = [0];
    let mut x_val: i8;
    let mut y_val: i8;
    let mut normalized_x_accel: f32;
    let mut normalized_y_accel: f32;
    let mut y_counter: usize = 0;
    let mut frame_count = 0;
    loop {
        Timer::after_millis(5).await;
        // if frame_count > 14 {
        //     embassy_rp::rom_data::reboot(0x0002, 1, 0x00, 0x01); // reboot to BOOTSEL
        // } else {
        //     frame_count += 1;
        // }
        i2c.write_read_async(addr, [0x29], &mut xh).await.unwrap(); // read accel data
        i2c.write_read_async(addr, [0x2B], &mut yh).await.unwrap(); // read accel data
        x_val = xh[0] as i8;
        y_val = yh[0] as i8;
        normalized_x_accel = (0.7 * accel_scale_max * (x_val as f32) / 128.0) / dt;
        normalized_y_accel = (0.7 * accel_scale_max * (y_val as f32) / 128.0) / dt;

        if y_val > 10 {
            y_counter += 1;
            if y_counter > 1000 {
                // embassy_rp::rom_data::reboot(0x0002, 1, 0x00, 0x01); // reboot to BOOTSEL
                Timer::after_millis(3000).await;
            }
        } else {
            if y_counter > 0 {
                y_counter -= 1;
            }
        }

        ACCEL_DATA_SIGNAL.signal([normalized_x_accel, -normalized_y_accel]);
    }
}

#[embassy_executor::task]
async fn simulation_update() {
    let mut scene = Scene::setupScene(500);
    ACCEL_DATA_SIGNAL.wait().await;
    let mut frame_count = 0;
    let mut miss_count = 0;
    let mut shake_count = 0;
    loop {
        if let First(accel_measurment) =
            select(ACCEL_DATA_SIGNAL.wait(), Timer::after_millis(10)).await
        {
            miss_count = 0;
            scene.set_gravity(accel_measurment);
            if accel_measurment[1] > 1000.0 || accel_measurment[1] < -1000.0 {
                shake_count += 100;
                if shake_count > 400 {
                    scene.particle_add(500, 500);
                }
            } else if shake_count > 0 {
                shake_count -= 1;
            }
        } else {
            miss_count += 1
        }
        if frame_count > 1000 {
            if frame_count % 10 == 0 {
                scene.particle_add(-1, 500);
            }
        }

        if !scene.is_paused() && scene.get_num_particles() > 0 && miss_count < 10 {
            frame_count += 1;
            scene.simulate();
            FRAME_DATA_SIGNAL.signal(scene.get_output());
        } else {
            // ACCEL_DATA_SIGNAL.wait().await;
        }
    }
}

fn setup_pio_2(
    pio: embassy_rp::Peri<'static, PIO0>,
    pin0: embassy_rp::Peri<'static, PIN_0>,
    pin1: embassy_rp::Peri<'static, PIN_1>,
    pin2: embassy_rp::Peri<'static, PIN_2>,
    pin3: embassy_rp::Peri<'static, PIN_3>,
    pin4: embassy_rp::Peri<'static, PIN_4>,
    pin5: embassy_rp::Peri<'static, PIN_5>,
    pin6: embassy_rp::Peri<'static, PIN_6>,
    pin7: embassy_rp::Peri<'static, PIN_7>,
    pin8: embassy_rp::Peri<'static, PIN_8>,
    pin9: embassy_rp::Peri<'static, PIN_9>,
    pin10: embassy_rp::Peri<'static, PIN_10>,
    pin11: embassy_rp::Peri<'static, PIN_11>,
    pin12: embassy_rp::Peri<'static, PIN_12>,
    pin13: embassy_rp::Peri<'static, PIN_13>,
    pin14: embassy_rp::Peri<'static, PIN_14>,
    pin15: embassy_rp::Peri<'static, PIN_15>,
    pin16: embassy_rp::Peri<'static, PIN_16>,
    pin17: embassy_rp::Peri<'static, PIN_17>,
    pin18: embassy_rp::Peri<'static, PIN_18>,
    pin19: embassy_rp::Peri<'static, PIN_19>,
    pin20: embassy_rp::Peri<'static, PIN_20>,
    pin21: embassy_rp::Peri<'static, PIN_21>,
) -> embassy_rp::pio::StateMachine<'static, PIO0, 0> {
    let Pio {
        mut common,
        sm0: mut sm,
        ..
    } = Pio::new(pio, Irqs); //contains the pio peripheral
    let prg = pio_asm!(
        //this program  is for moving data from TX to RX
        "out pins, 32", //change the state of pins, 0 is off, 1 is on.
        "out pindirs, 32",
        "nop [4]",
        "out pindirs, 32",
    );
    let mut cfg = Config::default(); // how does this know it's a PIO config?

    let p0 = common.make_pio_pin(pin0);
    let p1 = common.make_pio_pin(pin1);
    let p2 = common.make_pio_pin(pin2);
    let p3 = common.make_pio_pin(pin3);
    let p4 = common.make_pio_pin(pin4);
    let p5 = common.make_pio_pin(pin5);
    let p6 = common.make_pio_pin(pin6);
    let p7 = common.make_pio_pin(pin7);
    let p8 = common.make_pio_pin(pin8);
    let p9 = common.make_pio_pin(pin9);
    let p10 = common.make_pio_pin(pin10);
    let p11 = common.make_pio_pin(pin11);
    let p12 = common.make_pio_pin(pin12);
    let p13 = common.make_pio_pin(pin13);
    let p14 = common.make_pio_pin(pin14);
    let p15 = common.make_pio_pin(pin15);
    let p16 = common.make_pio_pin(pin16);
    let p17 = common.make_pio_pin(pin17);
    let p18 = common.make_pio_pin(pin18);
    let p19 = common.make_pio_pin(pin19);
    let p20 = common.make_pio_pin(pin20);
    let p21 = common.make_pio_pin(pin21);

    cfg.set_out_pins(&[
        &p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8, &p9, &p10, &p11, &p12, &p13, &p14, &p15, &p16,
        &p17, &p18, &p19, &p20, &p21,
    ]);
    cfg.use_program(&common.load_program(&prg.program), &[]); // load the program instructions
    cfg.clock_divider = (U56F8!(125_000_000) / 200_000).to_fixed(); //set the speed, I thing its 12,500Hz
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

    sm // this is where the state machine interface is
}

#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! {
    // embassy_rp::rom_data::reboot(0x0002, 1, 0x00, 0x01); // reboot to BOOTSEL
    loop {}
}
