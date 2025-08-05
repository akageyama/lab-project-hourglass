/*
  bounce_balls03_N_balls.pde
 
 - Show balls bouncing in the window. 
 - Click to toggle states of run/stop.
 */


boolean runningStateToggle = true;

float simulationTime = 0.0; 

int NBALLS = 10; 
Ball[] balls = new Ball[NBALLS];

class Ball {
  float x; // 粒子位置 x座標  // Ball's position x
  float y; // 粒子位置 y座標  //        position y
  float vx, vy; // 粒子速度 x成分とy成分 // Ball's velocity, x and y components.
  color col;  // 粒子の色 // color.

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
    stroke( 0 );
    fill( col );
    ellipse( x, y, 10, 10 );  // 円を描く // place a circle
  }

  void move(float dt) {    
    x += vx*dt;  // 粒子の移動 x方向 // shift the ball position in x.
    y += vy*dt;  // 粒子の移動 y方向 //  ... in y.
    if ( x >= width-1 && vx > 0 ) { // 粒子が右の壁に当たった // The ball hits the right wall.
           // ウィンドウはwidth個のピクセル分の幅をもつ。0からwidth-1までがxの取りうる範囲。
      vx = -vx; // 速度反転 // Reverse the direction.
    }
    if ( x <= 0 && vx < 0 ) {  // 粒子が左の壁に当たった// Hit the left wall.
      vx = -vx;
    }
    if ( y >= height-1 && vy > 0 ) { // 上の壁に衝突 // Hit the floor. 
      vy = -vy; // 速度反転 // Reverse the direction.
    }
    if ( y <= 0 && vy < 0 ) {  // 下の壁に衝突 // Hit the ceiling.
      vy = -vy;  // 速度反転  
    }
  }
}


void draw() {
  float dt = 0.1; // delta t (time increment).
  
  background( 255 );  

  if ( runningStateToggle ) {
    for ( int i=0; i<NBALLS; i++ ) {
      balls[i].move( dt );
    }  
    simulationTime += dt;  
    println( " t = " + simulationTime );
  } 

  for ( int i=0; i<NBALLS; i++ ) {
    balls[i].show();
  }
}


void setup() {
  size( 700, 700 );
  
  float x0 = width*0.4;
  float y0 = height*0.45;
  float velocity = 10.0;
  float radius = width*0.2;
  
  for ( int i=0; i<NBALLS; i++ ) {
    float theta = (TWO_PI/NBALLS)*i;
    // 初期位置 // initial position
    float x = x0 + radius*cos( theta );
    float y = y0 + radius*sin( theta );
    // 初期速度 // initial velocity
    float vx = velocity*cos( theta );
    float vy = velocity*sin( theta );
    // color
    int r = int( random( 100, 255 ) );
    int g = int( random( 100, 255 ) );
    int b = int( random( 100, 255 ) );
    balls[i] = new Ball( x, y, vx, vy, r, g, b ); 
  }
}


void mousePressed() {
  runningStateToggle = !runningStateToggle;
  if ( !runningStateToggle ) { // stopped.
    println( " Stopped." );
  }
}
