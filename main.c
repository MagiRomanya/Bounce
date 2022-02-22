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
  particle.r = 20.0f;
  particle.lastCollision = -1;
  particle.id = *counter;
  (*counter)++;
  return particle;
}

#define ScreenHeight 600
#define ScreenWidth 800


void drawParticles(Color color, Particle particles[], int n){
  for (int i =0; i<n; i++){
    DrawCircle(particles[i].x, particles[i].y, particles[i].r, color);
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

void resolveCollision(Particle* p1, Particle* p2){
  float vx1, vy1;
  float vx2, vy2;
  float dot = (p1->vx-p2->vx)*(p1->x-p2->x) + (p1->vy-p2->vy)*(p1->y-p2->y);
  float distance = (p1->r+p2->r)*(p1->r+p2->r);
  float tmass = p1->m+p2->m;
  float frac = dot/(distance*tmass);

  // Calculate the new velocities
  vx1 = p1->vx - 2.0f*p2->m*frac*(p1->x-p2->x);
  vy1 = p1->vy - 2.0f*p2->m*frac*(p1->y-p2->y);

  vx2 = p2->vx - 2.0f*p1->m*frac*(p2->x-p1->x);
  vy2 = p2->vy - 2.0f*p1->m*frac*(p2->y-p1->y);

  // Assign the new velocities
  p1->vx = vx1;
  p1->vy = vy1;
  p2->vx = vx2*1.0f;
  p2->vy = vy2*1.0f;
  printf("Colliding!\n");

}

float collisionDistance(float t, Particle p1, Particle p2){
  float output = 0;
  float x1t, x2t, y1t, y2t;
  // Calculate the position as a function of time
  x1t = p1.vx*t;
  y1t = p1.vy*t;
  x2t = p2.vx*t;
  y2t = p2.vy*t;
  output = (x1t-x2t)*(x1t-x2t) + (y1t-y2t)*(y1t-y2t);
  return output - (p1.r+p2.r)*(p1.r+p2.r);
}

float collisionWallTime(Particle p1){
  float epsilon = 0.001f; // Precission
  float xn[3] = {0.0f, 0.01f, 0.0f};
  bool found = false;
  float denominator;
  for (int i=0; i<1000; i++) {
    denominator = (xn[1]-xn[0]);
    xn[2] = xn[1] - xn[1]*(xn[1]-xn[0])/denominator;
    if (denominator < epsilon) printf("Denominator = %f\n",denominator);
    // Cycles through the array
    xn[0] = xn[1];
    xn[1] = xn[2];
    xn[2] = 0.0f;
    if (xn[1]<epsilon) {
      found = true;
      break;
    }
  }
  if (!found) printf("WARNING: Root not found! %f   ",xn[1]);
  return xn[1];
}

float collisionTime(Particle p1, Particle p2){
  float epsilon = 0.001f; // Precission
  float xn[3] = {0.0f, 0.01f, 0.0f};
  bool found = false;
  float denominator;
  for (int i=0; i<1000; i++) {
    denominator = (collisionDistance(xn[1],p1,p2)-collisionDistance(xn[0],p1,p2));
    xn[2] = xn[1] - collisionDistance(xn[1],p1,p2)*(xn[1]-xn[0])/denominator;
    if (denominator < epsilon) printf("Denominator = %f\n",denominator);
    // Cycles through the array
    xn[0] = xn[1];
    xn[1] = xn[2];
    xn[2] = 0.0f;
    if (collisionDistance(xn[1],p1,p2)<epsilon) {
      found = true;
      break;
    }
  }
  if (!found) printf("WARNING: Root not found! %f   ",xn[1]);
  return xn[1];
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

  for (int i=0; i<n; i++){
    // Collisions with screen edges
    if ((particles[i].x - particles[i].r<0) | (particles[i].x+particles[i].r>ScreenWidth)){
      particles[i].vx *= -1.0f;
    }
    if ((particles[i].y - particles[i].r<0) | (particles[i].y+particles[i].r>ScreenHeight)){
      particles[i].vy *= -1.0f;
    }
    // Collisions between other particles

    for (int j=0; j<i; j++){
      float futureX[4] ={
        particles[i].x+ particles[i].vx,
        particles[i].y+ particles[i].vy,
        particles[j].x+ particles[j].vx,
        particles[j].y+ particles[j].vy,
      };
      if (fabs(futureX[0]-futureX[2])<particles[i].r+particles[j].r) {
        if (fabs(futureX[1]-futureX[3])<particles[i].r+particles[j].r) {
          float dx = futureX[0]-futureX[2];
          float dy = futureX[1]-futureX[3];
          float dist = (dx)*(dx) + (dy)*(dy);
          float colDist = (particles[i].r+particles[j].r)*(particles[i].r+particles[j].r);
          if (dist < colDist) {
            // Buscar el punt on xoquen i calcular a on han de ser al següent instant de temps
            float time = collisionTime(particles[i], particles[j]);
            updateParticle(&particles[i],time);
            updateParticle(&particles[j],time);
            resolveCollision(&particles[i],&particles[j]);
            updateParticle(&particles[i],1.0f-time);
            updateParticle(&particles[j],1.0f-time);
          }
          else {
            updateParticle(&particles[i],1.0);
            updateParticle(&particles[j],1.0);
          }
        }
        else {
          updateParticle(&particles[i],1.0);
          updateParticle(&particles[j],1.0);
        }
      }
      else {
        updateParticle(&particles[i],1.0);
        updateParticle(&particles[j],1.0);
      }
      if (fabs(particles[i].x-particles[j].x)<particles[i].r+particles[j].r) {
        if (fabs(particles[i].y-particles[j].y)<particles[i].r+particles[j].r) {
          float dx = particles[i].x-particles[j].x;
          float dy = particles[i].y-particles[j].y;
          float dist = (dx)*(dx) + (dy)*(dy);
          float colDist = (particles[i].r+particles[j].r)*(particles[i].r+particles[j].r);
          if (dist < colDist) {
            if (!(particles[i].id == particles[j].lastCollision)) {
//              resolveCollision(&particles[i], &particles[j]);
              particles[i].lastCollision = particles[j].id;
              particles[j].lastCollision = particles[i].id;
            }
          else {
            particles[i].lastCollision = -1;
            particles[j].lastCollision = -1;
          }
          }
        }
      }
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
  particles[0] = genParticle(50.0f,50.0f,1.0f,1.0f, &Nparticles);
  particles[1] = genParticle(50.0f,150.0f,-1.5f,0.0f, &Nparticles);
  particles[2] = genParticle(150.0f,150.0f,1.0f,-1.0f, &Nparticles);
  particles[3] = genParticle(150.0f,350.0f,2.0f,0.5f, &Nparticles);
  bool pause = false;
  while(!WindowShouldClose()){
    if (IsKeyPressed(KEY_SPACE)) pause = !pause;
    if (!pause){
    // Update
    manageCollisions(particles, Nparticles);
    manageCollisions(particles, Nparticles);
    //updatePhysics(particles, Nparticles);
    }
    // Drawing
    BeginDrawing();
      ClearBackground(RAYWHITE);
      drawParticles(RED, particles, Nparticles);
      drawPause(pause);
    EndDrawing();
  }
  CloseWindow(); // Closes the window
  return 0;
}
