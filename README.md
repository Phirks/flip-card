This repo contains all files related to the flip-card project, which is a business card that runs a fluid-implicit-particle(FLIP) simulation.

![alt text](media/1000003136.gif)

The fluid simulation logic is contained in a standalone crate, which is in the "fluid_sim_crate" folder

a WASM simulator is also provided in the "sim_display" folder, which is what I use to debug issues in the simulation.

The implementation for the fluid simulation on the rp2350 is in the "flip-card_firmware" file

further details can be found in each folder's README files
