#include "raylib.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Particle data structure
typedef struct{
  int id; // Particle identifier
  float x; // Horizontal coordinate
  float y; // Vertical coordinate
  float vx; // Horizontal velocity component
  float vy; // Vertical velocity component
  float r; // Particle's radious
  float m; // Particle's mass
  int lastCollision; // Id of the last particle it collided with
}Particle;


Particle genParticle(float x, float y, float vx, float vy, int* counter){
  Particle particle;
  particle.x = x;
  particle.y = y;
  particle.vx = vx;
  particle.vy = vy;
  particle.m = 1.0f;
  particle.r = 20.0;
  particle.lastCollision = -1;
  particle.id = *counter;
  (*counter)++;
  return particle;
}

#define ScreenHeight 600
#define ScreenWidth 800
#define DeltaT 1.0f


void drawParticles(Particle particles[], int n){
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
  if (n>n_colors) {
    printf("WARNING: Too much particles for the colors");
  }
  for (int i =0; i<n; i++){
    DrawCircle(particles[i].x, particles[i].y, particles[i].r, colors[i]);
    /* DrawCircle(particles[i].x + particles[i].vx*DeltaT, */
    /*            particles[i].y + particles[i].vy*DeltaT, */
    /*            particles[i].r, colors[i]); */
    // Draw velocities
    //DrawLine(particles[i].x, particles[i].y, (particles[i].vx*10+particles[i].x), (particles[i].vy*10+particles[i].y), BLACK);
  }
}

void drawPause(bool pause){
  Vector2 coordinates = (Vector2){ScreenWidth-50, 50};
  if (pause) {
    Vector2 size = (Vector2){-3, -20};
    Color pauseColor = LIGHTGRAY;
    DrawRectangleV(coordinates,size, pauseColor);
    DrawRectangleV((Vector2){coordinates.x+10, coordinates.y},size, pauseColor);
  }
}

void drawStatistics(Particle particles[], int n){
  float tmomentum = 0.0f;
  float tenergy = 0.0f;
  float speed;
  for (int i=0;i<n;i++){
    speed = sqrtf(particles[i].vx*particles[i].vx+particles[i].vy*particles[i].vy);
    tmomentum += particles[i].m*speed;
    tenergy += particles[i].m*speed*speed/2.0f;
  }
  DrawText(TextFormat("Momentum = %f",tmomentum), 40, 50, 20, LIGHTGRAY);
  DrawText(TextFormat("Energy = %f",tenergy), 40, 80, 20, LIGHTGRAY);
}
void resolveCollision(Particle* p1, Particle* p2){
  float vx1, vy1;
  float vx2, vy2;
  float dot = (p1->vx-p2->vx)*(p1->x-p2->x) + (p1->vy-p2->vy)*(p1->y-p2->y);
  //float distance = (p1->r+p2->r)*(p1->r+p2->r);
  float distance = (p1->x-p2->x)*(p1->x-p2->x) + (p1->y-p2->y)*(p1->y-p2->y);
  float tmass = p1->m+p2->m;
  float frac = dot/(distance*tmass);

  // Calculate the new velocities
  vx1 = p1->vx - 2.0f*p2->m*frac*(p1->x-p2->x);
  vy1 = p1->vy - 2.0f*p2->m*frac*(p1->y-p2->y);

  vx2 = p2->vx - 2.0f*p1->m*frac*(p2->x-p1->x);
  vy2 = p2->vy - 2.0f*p1->m*frac*(p2->y-p1->y);

  //printf("Colliding!   D=%f, Delta=%f\n",distance,distance - (p1->r+p2->r)*(p1->r+p2->r));
  // Assign the new velocities
  p1->vx = vx1;
  p1->vy = vy1;
  p2->vx = vx2;
  p2->vy = vy2;
}

float collisionDistance(float t, Particle p1, Particle p2){
  float output = 0;
  float x1t, x2t, y1t, y2t;
  // Calculate the position as a function of time
  x1t = p1.x + p1.vx*t;
  y1t = p1.y + p1.vy*t;
  x2t = p2.x + p2.vx*t;
  y2t = p2.y + p2.vy*t;
  output = (x1t-x2t)*(x1t-x2t) + (y1t-y2t)*(y1t-y2t);
  return output - (p1.r+p2.r)*(p1.r+p2.r);
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

float zeroSecant(float func(float), float guess, float guess2, float epsilon){
  float xn[3] = {guess, guess2, 0.0f};
  bool found = false;
  for (int i=0; i<1000; i++) {
    xn[2] = xn[1] - func(xn[1])*(xn[1]-xn[0])/(func(xn[1])-func(xn[0]));
    xn[0] = xn[1];
    xn[1] = xn[2];
    xn[2] = 0.0f;
    if (func(xn[1])<epsilon) {
      found = true;
      break;
    }
  }
  if (!found) printf("WARNING: Root not found!");
  return xn[1];
}

void updateParticle(Particle* particle, float t){
  particle->x += particle->vx*t;
  particle->y += particle->vy*t;
}

void manageCollisions(Particle particles[], int n){
  bool collided[n];
  float nextX[2];
  for (int i=0; i<n; i++){
    nextX[0] = particles[i].x+ particles[i].vx*DeltaT;
    nextX[1] = particles[i].y+ particles[i].vy*DeltaT;
    // Collisions with screen edges
    if ((nextX[0] - particles[i].r<0) | (nextX[0]+particles[i].r>ScreenWidth)){
      collided[i] = true;
      particles[i].vx *= -1.0f;
    }
    if ((nextX[1] - particles[i].r<0) | (nextX[1]+particles[i].r>ScreenHeight)){
      collided[i] = true;
      particles[i].vy *= -1.0f;
    }
    // Collisions between other particles

    for (int j=0; j<i; j++){
      float futureX[4] ={
        particles[i].x+ particles[i].vx*DeltaT,
        particles[i].y+ particles[i].vy*DeltaT,
        particles[j].x+ particles[j].vx*DeltaT,
        particles[j].y+ particles[j].vy*DeltaT,
      };
      if (fabs(futureX[0]-futureX[2])<particles[i].r+particles[j].r) {
        if (fabs(futureX[1]-futureX[3])<particles[i].r+particles[j].r) {
          float dx = futureX[0]-futureX[2];
          float dy = futureX[1]-futureX[3];
          float dist = (dx)*(dx) + (dy)*(dy);
          float colDist = (particles[i].r+particles[j].r)*(particles[i].r+particles[j].r);
          float dxp = particles[i].x-particles[j].x;
          float dyp = particles[i].y-particles[j].y;
          if ((dist < colDist) && !(dxp*dxp + dyp*dyp < colDist)) {
            // Buscar el punt on xoquen i calcular a on han de ser al següent instant de temps
            float time = collisionTime(particles[i], particles[j]);
            collided[i] = true;
            collided[j] = true;
            updateParticle(&particles[i],time*DeltaT);
            updateParticle(&particles[j],time*DeltaT);
            resolveCollision(&particles[i],&particles[j]);
            updateParticle(&particles[i],DeltaT*(1.0f-time));
            updateParticle(&particles[j],DeltaT*(1.0f-time));
          }
        }
      }
    }
    if (!collided[i]) {
      updateParticle(&particles[i], DeltaT);
    }
  }
}

void updatePhysics(Particle particles[], int n){
  for (int i=0; i<n; i++){
    particles[i].x += particles[i].vx;
    particles[i].y += particles[i].vy;
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
    if (IsKeyPressed(KEY_SPACE)) pause = !pause;
    if (!pause){
    // Update
    manageCollisions(particles, Nparticles);
    }
    else if (IsKeyPressed(KEY_RIGHT)) {
      manageCollisions(particles, Nparticles);
    }
    // Drawing
    BeginDrawing();
      ClearBackground(RAYWHITE);
      drawParticles(particles, Nparticles);
      drawPause(pause);
      drawStatistics(particles,Nparticles);
    EndDrawing();
  }
  CloseWindow(); // Closes the window
  return 0;
}
