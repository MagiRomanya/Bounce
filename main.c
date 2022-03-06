#include "raylib.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float InvSqrt(float x)
{
        float xhalf = 0.5f * x;
        int i = *(int*)&x;            // store floating-point bits in integer
        i = 0x5f3759df - (i >> 1);    // initial guess for Newton's method
        x = *(float*)&i;              // convert new bits into float
        x = x*(1.5f - xhalf*x*x);     // One round of Newton's method
        return x;
}
// Particle data structure
typedef struct{
  int id; // Particle identifier
  float x; // Horizontal coordinate
  float y; // Vertical coordinate
  float vx; // Horizontal velocity component
  float vy; // Vertical velocity component
  float x1, y1, vx1, vy1; // Future coordinates and velocities
  float r; // Particle's radious
  float m; // Particle's mass
}Particle;


Particle genParticle(float x, float y, float vx, float vy, int* counter){
  /* Initializes a particle struct */
  Particle particle;
  particle.x = x;
  particle.y = y;
  particle.vx = vx;
  particle.vy = vy;
  // We calculate the future positions and velocities later
  particle.x1 = x;
  particle.y1 = y;
  particle.vx1 = vx;
  particle.vy1 = vy;
  particle.m = 1.0f;
  particle.r = 20.0;
  particle.id = *counter;
  (*counter)++;
  return particle;
}

#define ScreenHeight 600
#define ScreenWidth 800
#define DeltaT 0.1f
// TODO: Fix the overlapping particles: (BUG: calculating new coldinates with old velocities)

void drawParticles(Particle particles[], int n){
  /* Draws the particles to the screen in diferent colors */
  // TODO: Make the colors cycle when Nparticles > Ncolors
  Color colors[]= {
    DARKPURPLE,
    PINK,
    RED,
    MAROON,
    GREEN,
    LIME,
    BLUE,
    DARKBLUE,
    SKYBLUE,
    PURPLE,
    VIOLET
  };
  int n_colors=sizeof(colors)/sizeof(RED);
  /* for (int i =0; i<n; i++){ */
  /*   DrawCircle(particles[i].x1, particles[i].y1, particles[i].r, BLACK); */
  /* } */
  for (int i =0; i<n; i++){
    DrawCircle(particles[i].x, particles[i].y, particles[i].r, colors[i%n_colors]);
    // Draw velocities
    //DrawLine(particles[i].x, particles[i].y, (particles[i].vx*10+particles[i].x), (particles[i].vy*10+particles[i].y), BLACK);
  }
}

void drawPause(bool pause){
  /* Dispays the pause symbol when the game is paused */
  Vector2 coordinates = (Vector2){ScreenWidth-50, 50};
  if (pause) {
    Vector2 size = (Vector2){-3, -20};
    Color pauseColor = LIGHTGRAY;
    DrawRectangleV(coordinates,size, pauseColor);
    DrawRectangleV((Vector2){coordinates.x+10, coordinates.y},size, pauseColor);
  }
}

void drawStatistics(Particle particles[], int n){
  /* Makes visible the energy and momentum counter on the screen */
  float tmomentum = 0.0f;
  float tenergy = 0.0f;
  float speed2;
  for (int i=0;i<n;i++){
    speed2 = (particles[i].vx*particles[i].vx+particles[i].vy*particles[i].vy);
    tmomentum += particles[i].m*sqrtf(speed2);
    tenergy += particles[i].m*speed2/2.0f;
  }
  DrawText(TextFormat("Momentum = %f",tmomentum), 40, 50, 20, LIGHTGRAY);
  DrawText(TextFormat("Energy = %f",tenergy), 40, 80, 20, LIGHTGRAY);
}
void resolveCollision(Particle* p1, Particle* p2){
  /* Handles the physics of a collision between 2 particles (elastic collisions). It assigns
   *the new velocities to the particles */
  float vx_new1, vy_new1;
  float vx_new2, vy_new2;
  float dot = (p1->vx1-p2->vx1)*(p1->x1-p2->x1) + (p1->vy1-p2->vy1)*(p1->y1-p2->y1);
  //float distance = (p1->r+p2->r)*(p1->r+p2->r);
  float distance = (p1->x1-p2->x1)*(p1->x1-p2->x1) + (p1->y1-p2->y1)*(p1->y1-p2->y1);
  float tmass = p1->m+p2->m;
  float frac = dot/(distance*tmass);

  // Calculate the new velocities
  vx_new1 = p1->vx - 2.0f*p2->m*frac*(p1->x1-p2->x1);
  vy_new1 = p1->vy - 2.0f*p2->m*frac*(p1->y1-p2->y1);

  vx_new2 = p2->vx - 2.0f*p1->m*frac*(p2->x1-p1->x1);
  vy_new2 = p2->vy - 2.0f*p1->m*frac*(p2->y1-p1->y1);

  //printf("Colliding!   D=%f, Delta=%f\n",distance,distance - (p1->r+p2->r)*(p1->r+p2->r));
  // Assign the new velocities
  p1->vx1 = vx_new1;
  p1->vy1 = vy_new1;
  p2->vx1 = vx_new2;
  p2->vy1 = vy_new2;
}

float collisionWallTime(Particle p1){
  float time;
  float dx, dy;
  dx = p1.x + p1.vx*DeltaT;
  dy = p1.y + p1.vy*DeltaT;
  if (dx + p1.r> ScreenWidth) {
    time = fabs(ScreenWidth-p1.x-p1.r)/p1.vx;
  }
  else if (dx - p1.r < 0) {
    time = (p1.x-p1.r)/p1.vx;
  }
  else if (dy + p1.r > ScreenHeight) {
    time = fabs(ScreenHeight-p1.y-p1.r)/p1.vy;
  }
  else if (dy - p1.r < 0) {
    time = (p1.y-p1.r)/p1.vy;
  }
  return time;
}

float collisionTime(Particle p1, Particle p2){
  /* Calculates when 2 particles would collide assuming they are not colliding and will have collided in the next frame */
  float dx, dy, dvx, dvy;
  float a, b, c; // 2nd degree equation coefficients
  float time;
  dx = p2.x - p1.x;
  dy = p2.y - p1.y;
  dvx = (p2.vx - p1.vx)*DeltaT;
  dvy = (p2.vy - p1.vy)*DeltaT;
  a = dvx*dvx + dvy*dvy;
  b = 2.0f*(dx*dvx + dy*dvy);
  c = dx*dx + dy*dy - (p1.r+p2.r)*(p1.r+p2.r);
  time = (-b-sqrt(b*b - 4.0f*a*c))/(2.0f*a);
  if (time < 0.0f || time > 1.0f) {
    printf("WARNING: Collision time outside frame");
  }
  return time;
}

void updateParticle(Particle* particle, float t){
  particle->x1 = particle->x + particle->vx*t;
  particle->y1 = particle->y + particle->vy*t;
}

void manageCollisions(Particle particles[], int n){
  /* Checks for collisions with the screen border and between particles */
  for (int i=0; i<n; i++){
    // Collisions with screen edges
    if ((particles[i].x1 - particles[i].r<0) | (particles[i].x1+particles[i].r>ScreenWidth)){
      particles[i].vx1 = -1.0f*particles[i].vx;
    }
    if ((particles[i].y1 - particles[i].r<0) | (particles[i].y1+particles[i].r>ScreenHeight)){
      particles[i].vy1 = -1.0f*particles[i].vy;
    }
    // Collisions between other particles

    for (int j=0; j<i; j++){
      if (fabs(particles[i].y1-particles[j].y1)<particles[i].r+particles[j].r) {
        if (fabs(particles[i].x1-particles[j].x1)<particles[i].r+particles[j].r) {
          float dx = particles[i].x1-particles[j].x1;
          float dy = particles[i].y1-particles[j].y1;
          float dist = (dx)*(dx) + (dy)*(dy);
          float colDist = (particles[i].r+particles[j].r)*(particles[i].r+particles[j].r);
          float dxp = particles[i].x-particles[j].x;
          float dyp = particles[i].y-particles[j].y;
          if ((dist < colDist) && !(dxp*dxp + dyp*dyp < colDist)) {
            // Buscar el punt on xoquen i calcular a on han de ser al següent instant de temps
            float time = collisionTime(particles[i], particles[j]);
            particles[i].x1 = particles[i].x + particles[i].vx*time*DeltaT;
            particles[i].y1 = particles[i].y + particles[i].vy*time*DeltaT;
            particles[j].x1 = particles[j].x + particles[j].vx*time*DeltaT;
            particles[j].y1 = particles[j].y + particles[j].vy*time*DeltaT;
            resolveCollision(&particles[i],&particles[j]);
            particles[i].x1 = particles[i].x1 + particles[i].vx1*(1.0f-time)*DeltaT;
            particles[i].y1 = particles[i].y1 + particles[i].vy1*(1.0f-time)*DeltaT;
            particles[j].x1 = particles[j].x1 + particles[j].vx1*(1.0f-time)*DeltaT;
            particles[j].y1 = particles[j].y1 + particles[j].vy1*(1.0f-time)*DeltaT;
          }
        }
      }
    }
  }
}

// Calculates a step of RK4
//
//	EXAMPLE: Harmonic oscilator
//
//	diferential eq  dv_z/dt = -kz     (where v_z = dz/dt)
//		        dz/dt = v_z
//
//	=> y = (v_z, z)	=> dy/dx = f(x,y) = (dv_z/dz, dz/dx) = ( -kz, v_z))
//	   x = t
// y[]       Dependent variables (position & velociteies)
// x         Independent variable
// size      Number of dependent variables
// func      Function which returns the value of the derivative
// h         Step size
// y1[]      Updated dependent variables for the next step
//
void rungekutta4(const float y[], float x, int size, int (*func)(float, const float*, float*, int), float h, float* y1)
{
	int i;
	float K1[size], K2[size], K3[size], K4[size], ycache[size];

	// Sets K1 = f(x,y)
	func(x,y,K1,size);

	// Sets K2 = f(x + h/2, y + h/2*K1)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h/2*K1[i];
	}
	func(x + h/2, ycache, K2, size);

	// Sets K3 = f(x + h/2, y + h/2*K2)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h/2*K2[i];
	}
	func(x + h/2, ycache, K3, size);

	// Sets K4 = f(x + h, y + h*K3)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h*K3[i];
	}
	func(x + h, ycache, K4, size);

	// Computes the final answer
	for (i=0; i<size; i++)
	{
		y1[i] = y[i] + h/6*(K1[i] + 2*(K2[i]+K3[i]) + K4[i]);
	}
}

int gravity(float x, const float y[], float f[], int size)
{
	int i, j, k;
	const int SIZE = size;
	const int SIZE2 = size/2; // SIZE/2
	const int N=SIZE/4; // We have N particles, DIM coordinates and DIM velocities => number of dimentions of y[] f[]
	const float G = 1.0e3f; // Gravitational constant multiplied by m
  const int DIM = 2;
	float sum2, dist;


	// Initialize f
	for (i=0; i<SIZE2; i++)
	{
		f[i] = y[SIZE2+i]; // dz/dx = v_z
		f[SIZE2+i] = 0;
	}
	for (i=0; i<SIZE2 -DIM; i+=DIM)
	{
		for (j=i+DIM; j<SIZE2; j+=DIM)
		{
			sum2 = 0;
			// CALCULATE THE DISTANCE
			for (k=0; k<DIM; k++)
			{
				sum2 +=(y[j+k]-y[i+k])*(y[j+k]-y[i+k]); // x**2 + y**2 + z**2
			}
			dist = InvSqrt(sum2);
			// COMPUTES THE ACCELERATION
			for (k=0; k<DIM; k++)
			{
				double f_module = G*(dist*dist*dist)*(y[j+k]-y[i+k]);
				f[SIZE2+i+k] += f_module;
				f[SIZE2+j+k] += -f_module;
			}
		}
	}
	return 0;
}

void calculateNextStep_rk4(Particle particles[], int n,float time, float dt){
  float y[n*4];
  float y1[n*4];
  // Creates the y[] output vector
  for (int i=0; i<n; i++){
    y[i*2] = particles[i].x;
    y[i*2+1] = particles[i].y;
  }
  for (int i=0; i<n; i++){
    y[n*2+i*2] = particles[i].vx;
    y[n*2+i*2+1] = particles[i].vy;
  }
  rungekutta4(y, time, 4*n, gravity, DeltaT, y1);
  // Assigns the y1[] to the future coordinates and velocities
  for (int i=0; i<n; i++){
    particles[i].x1 = y1[i*2];
    particles[i].y1 = y1[i*2+1];
  }
  for (int i=0; i<n; i++){
    particles[i].vx1 = y1[n*2+i*2];
    particles[i].vy1 = y1[n*2+i*2+1];
  }
}

void calculateNextStep(Particle particles[], int n, float dt){
  /* This function must update the future positions and velocites according to the interactions and physical state */
  for (int i=0; i<n; i++){
    particles[i].x1 = particles[i].x + particles[i].vx*dt;
    particles[i].y1 =  particles[i].y + particles[i].vy*dt;
    particles[i].vx1 = particles[i].vx;
    particles[i].vy1 =  particles[i].vy;
  }
}

void updateParticles(Particle particles[], int n){
  /* Sets the future coordinates and velocities to the present ones*/
  for (int i = 0; i<n; i++){
    particles[i].x = particles[i].x1;
    particles[i].y = particles[i].y1;
    particles[i].vx = particles[i].vx1;
    particles[i].vy = particles[i].vy1;
  }
}

int main(void){
  InitWindow(ScreenWidth, ScreenHeight, "Bounce"); // Creates a window
  SetTargetFPS(60);

  // Initial particle distribution
  Particle particles[30];
  int Nparticles = 0;
  //                          x     y     vx   vy   r     m
  particles[0] = genParticle(50.0f,50.0f,1.5f,3.5f, &Nparticles);
  particles[1] = genParticle(50.0f,150.0f,-0.75f,3.0f, &Nparticles);
  particles[2] = genParticle(150.0f,150.0f,1.5f,-3.5f, &Nparticles);
  particles[3] = genParticle(150.0f,350.0f,3.5f,0.25f, &Nparticles);
  particles[4] = genParticle(50.0f,550.0f,-0.75f,3.0f, &Nparticles);
  particles[5] = genParticle(150.0f,550.0f,1.5f,-3.5f, &Nparticles);
  particles[6] = genParticle(150.0f,450.0f,3.5f,0.25f, &Nparticles);
  bool pause = false;

  while(!WindowShouldClose()){
    // Pauses the game
    if (IsKeyPressed(KEY_SPACE)) pause = !pause;
    if (!pause){
    // Update
      calculateNextStep_rk4(particles, Nparticles,0.0f, DeltaT);
      manageCollisions(particles, Nparticles);
    }
    else if (IsKeyPressed(KEY_RIGHT)) {
      calculateNextStep_rk4(particles, Nparticles,0.0f, DeltaT);
      manageCollisions(particles, Nparticles);
    }
    // Drawing
    BeginDrawing();
      ClearBackground(RAYWHITE);
      drawParticles(particles, Nparticles);
      drawPause(pause);
      drawStatistics(particles,Nparticles);
    EndDrawing();
    if (!pause){
    updateParticles(particles,Nparticles);
    }
    else if (IsKeyPressed(KEY_RIGHT)) {
    updateParticles(particles,Nparticles);
    }
  }
  CloseWindow(); // Closes the window
  return 0;
}
