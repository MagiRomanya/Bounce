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
      if (fabs(particles[i].x-particles[j].x)<particles[i].r+particles[j].r) {
        if (fabs(particles[i].y-particles[j].y)<particles[i].r+particles[j].r) {
          float dx = particles[i].x-particles[j].x;
          float dy = particles[i].y-particles[j].y;
          float dist = (dx)*(dx) + (dy)*(dy);
          float colDist = (particles[i].r+particles[j].r)*(particles[i].r+particles[j].r);
          if (dist < colDist) {
            if (!(particles[i].id == particles[j].lastCollision)) {
              resolveCollision(&particles[i], &particles[j]);
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
  particles[0] = genParticle(50.0f,50.0f,2.0f,2.0f, &Nparticles);
  particles[1] = genParticle(50.0f,150.0f,-4.0f,0.0f, &Nparticles);
  particles[2] = genParticle(150.0f,150.0f,4.0f,1.0f, &Nparticles);
  particles[3] = genParticle(150.0f,350.0f,4.0f,1.0f, &Nparticles);
  bool pause = false;
  while(!WindowShouldClose()){
    if (IsKeyPressed(KEY_SPACE)) pause = !pause;
    if (!pause){
    // Update
    manageCollisions(particles, Nparticles);
    updatePhysics(particles, Nparticles);
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
