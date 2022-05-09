!~!~For steadyFlow, we define array up, vp, Tp to check  convergence;
!~!~for unsteadyFlow, we did not define array up, vp, Tp to save memory.
#define steadyFlow    
! #define unsteadyFlow

!~!~Uncomment below to simulate mass particles
!~#define pointParticle

!!!~~velocity B.C.~~
#define HorizontalWallsNoslip
#define VerticalWallsNoslip
!!!~#define VerticalWallsPeriodicalU
!!!~#define HorizontalWallsFreeslip
!!!~#define VerticalWallsFreeslip
!!!~~velocity B.C.~~

!!!!~~temperature B.C. (for Rayleigh Benard Cell)~~
! #define RayleighBenardCell
! #define HorizontalWallsConstT
! #define VerticalWallsAdiabatic
!~#define VerticalWallsPeriodicalT
!!!!~~temperature B.C.~~

!!!~~temperature B.C. (for Side Heated Cell)~~
#define SideHeatedCell
#define HorizontalWallsAdiabatic
#define VerticalWallsConstT
!!!~~temperature B.C.~~