/*
  freefall_balls_rungekutta.pde
 
 - Show balls bouncing in the window. 
 - Click to toggle states of run/stop.
 */


boolean runningStateToggle = true;
float simulationTime = 0.0;
int simulationStep = 0;

int NBALLS = 3; 
float BALL_MASS = 1.0; 
float GRAVITY_ACCELERATION = 9.8;


float energy0;

Ball[] balls = new Ball[NBALLS];


class Ball {
  float x, y;      // Ball's position x & y.
  float vx, vy; // Ball's velocity, x & y components.
  color col;    // Ball's color.

  Ball( float xInit, float yInit, 
        float vxInit, float vyInit, 
        int red, int green, int blue ) {
    x = xInit;
    y = yInit;
    vx = vxInit;
    vy = vyInit;
    col = color( red, green, blue );
  }
  
  void show() {
    stroke( 0 );  // black circumference ring
    fill( col );
    ellipse( x, y, 10, 10 );  // place a circle
  }

  float getEnergy() {
    float vsq = pow( vx, 2 ) + pow( vy, 2 );   // velocity squared
    float kinetic_energy = 0.5*BALL_MASS*vsq;
    float potential_energy = BALL_MASS * GRAVITY_ACCELERATION * ( height - y );
    return kinetic_energy + potential_energy;
  }

  void reflect( float[] posvel ) {
    if ( ( posvel[0] > width ) && posvel[2] > 0 ) {// Hit the right wall.
      posvel[2] = -posvel[2]; // Reverse the direction.
    }
    if ( ( posvel[0] < 0 ) && posvel[2] < 0 ) {  // Hit the left wall.
      posvel[2] = -posvel[2];
    }
    if ( ( posvel[1] > height ) && posvel[3] > 0 ) {// Hit the floor.
      posvel[3] = -posvel[3]; // Reverse the direction.
    }
    if ( ( posvel[1] < 0 ) && posvel[3] < 0 ) {  // Hit the ceiling.
      posvel[3] = -posvel[3];
    }
  }

  void move( float dt ) {    
    float[] posvel = new float[4];
    posvel[0] = x;
    posvel[1] = y;
    posvel[2] = vx;
    posvel[3] = vy;

    //integrator__Euler( posvel, dt );
    integrator__RungeKutta4( posvel, dt ); 

    reflect( posvel );

    x = posvel[0];
    y = posvel[1];
    vx = posvel[2];
    vy = posvel[3];
  }
} 


float totalEnergy() {
  float sum = 0.0;
  for( int i=0; i<NBALLS; i++ ) {
    sum += balls[i].getEnergy();
  }
  return sum;
}


void simulation_step( float dt )
{
  float energy, error;
  
  if ( runningStateToggle ) {
    for ( int i=0; i<NBALLS; i++ ) {
      balls[i].move( dt );
    }  
    
    if ( simulationStep % 100 == 0 ) {
      energy = totalEnergy();
      error = abs( ( energy-energy0 ) / energy0 ); 
      println( " t=" + simulationTime 
             + " step=" + simulationStep
             + " energy=" + energy
             + " error="  + error );
    }
    simulationTime += dt;
    simulationStep++;
  }  
}


void draw() {
  background( 255 );  

  float dt = 0.01;   // delta t (time increment).

  for ( int s=0; s<10; s++ ) {
    simulation_step( dt );
  }

  for ( int i=0; i<NBALLS; i++ ) {
    balls[i].show();
  }
}


void setup() {
  size( 500, 400 );
  float maxVelocity = 50.0;
  float minVelocity =  0.0;
  
  for ( int i=0; i<NBALLS; i++ ) {
    // position
    float x0 = random( width );
    float y0 = random( height );
    // velocity
    float vel = random( minVelocity, maxVelocity );
    float angle = random( 0, TWO_PI );
    float vx0 = vel*cos( angle );
    float vy0 = vel*sin( angle );
    // color
    int r = int( random( 100, 255 ) );
    int g = int( random( 100, 255 ) );
    int b = int( random( 100, 255 ) );
    balls[i] = new Ball( x0, y0, vx0, vy0, r, g, b ); 
  }
  
  energy0 = totalEnergy();  // Total energy at t=0. 
}


void mousePressed() {
  runningStateToggle = !runningStateToggle;
  if ( !runningStateToggle ) { // stopped.
    // saveFrame();
    // println( " Stopped. Frame saved." );
  }
}

  
void equation_of_motion( float[] posvel, float dt, float[] delta ) 
{
  // The equation of motion of a particle is 
  //        dx/dt = vx,
  //        dy/dt = vy,
  //       dvx/dt = Force_x / mass,
  //       dvy/dt = Force_y / mass.
  // We combine the independent variables into an array as
  //        x --> posvel[0],
  //        y --> posvel[1],
  //       vx --> posvel[2],
  //       vy --> posvel[3].
  // The equation of motion is then written in generalized form:
  //       d(posvel[i])/dt = RHS[i]  (for i=0...3),
  // where RHS stands for right-hand-side of the equations.
  // Defining a new array delta as
  //       delta[i] = d(posvel[i]),
  // The equation of motion becomes
  //       delta[i] = RHS[i]*dt   (for i=0...3).
  //       
  delta[0] = posvel[2]*dt; // dx = vx*dt;
  delta[1] = posvel[3]*dt; // dy = vy*dt;
  delta[2] = 0;            // dvx = 0
  delta[3] = GRAVITY_ACCELERATION * dt; 
                           // dvy = G*dt
}

  
void integrator__Euler( float[] posvel, float dt ) 
{    
  float[] delta = new float[4];
  equation_of_motion( posvel, dt, delta );

  for ( int i=0; i<4; i++ ) {
    posvel[i] += delta[i];
  }
}


void integrator__RungeKutta4( float[] posvel, float dt ) 
{
  float[] delta1 = new float[4];
  float[] delta2 = new float[4];
  float[] delta3 = new float[4];
  float[] delta4 = new float[4];
  float[]   work = new float[4];
  final float ONE_SIXTH = 1.0/6.0;
  
  // Step 1
  equation_of_motion( posvel, dt, delta1 );
  runge_kutta4_advance( posvel, delta1, 0.5, work );
  
  // Step2
  equation_of_motion( work, dt, delta2 );
  runge_kutta4_advance( posvel, delta2, 0.5, work );
  
  // Step3
  equation_of_motion( work, dt, delta3 );
  runge_kutta4_advance( posvel, delta3, 1.0, work );

  // Step4
  equation_of_motion( work, dt, delta4 );

  // The result
  for ( int i=0; i<4; i++ ) {
    posvel[i] += (   delta1[i]
                 + 2*delta2[i]
                 + 2*delta3[i]
                 +   delta4[i] ) * ONE_SIXTH;
  }
}


void runge_kutta4_advance( float[] prev, float[] delta, 
                           float factor, float[] work ) 
{
  for ( int i=0; i<4; i++ ) {
    work[i] = prev[i] + factor*delta[i];
  }
}
