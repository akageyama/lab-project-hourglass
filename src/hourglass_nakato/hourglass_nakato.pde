/*

  hourglass_template.pde

         
         +---------+
         |         |
         |         |
         |         |
          \   .   /
           \  .  /
           /  .  \
          /   .   \  
         |    .    |
         |    .    |
         |         |
         +---------+
         
      
 【概要】
    * 砂時計 (hourglass) の1次元シミュレーションを実現するためのサンプル
    * ここでは砂粒 (sand grain) の自由落下と床面との衝突を解く
    * 鉛直方向(y方向)の1次元運動を解く
    * 複数（NSGIP個）の質点がy方向に1列に並ぶ
    * これを砂柱 (sand pillar）と呼ぶ
    * NSGIP=1とすると一つの質点の自由落下問題を解くシミュレーションになる
    * 砂柱の中の砂粒の順番が入れ替わる（ある砂粒が他を追い越す）ことはない
    * 砂柱の一番下の砂粒が落下し、床面と衝突すると床面から反発力を受ける
    * 衝突時にエネルギーの散逸があれば一番下の砂粒は最終的には床面上で静止する
    * 砂柱全体が静止状態のとき、砂柱の一番下の砂粒を床面が支える力が砂柱を支える垂直抗力
    * 床面が砂津橋らを支える垂直抗力から砂柱の重さを計算する
    * 実際には垂直抗力は時間的に激しく変化するので時間平均をとる
    * 砂柱の中の砂粒同士の相互作用は仮想的なバネを想定している
    
    
         |         o o
         |        o   o    <---  砂粒 1
         ↓         o o
         　         |
         重         |      <---  仮想バネ   
         力         |
         加        o o
         速       o   o    <---  砂粒 0 
         度        o o
                  
              -----------  <---  床面

 
 【計算手法】
    * 時間積分は4次のルンゲ・クッタ法を採用している
    * draw()関数の中からrungeKutta4()関数を呼ぶ
    * rungeKutta4()とrungeKutta4Advance()は別ファイルにある
    * Processingのエディタ画面では別のファイルは別のタブで表示される
    * rungeKutta4()とrungeKutta4Advance()を変更する必要はおそらくない
    * 砂粒と床面との相互作用にはバネ・ダンパーモデルを用いる


 【砂時計の1次元モデルへの改訂方針】
    * 砂粒間に想定している仮想バネを外し、砂粒同士の弾性衝突を導入する
    * 砂粒の運動を2次元または3次元的に追跡する必要はない
    * 床面をもう一つ、砂時計の真ん中の高さ (y=0) に設定する
    * 追加したy=0の床面に砂が落ちる穴（オリフィス）を設定する
    * オリフィスから定期的に砂粒が落ちるようにする
    * 2枚の床面がそれぞれ砂柱を支えるその垂直抗力の和を計算する
    * 上の床面から砂が落ちている最中にその垂直抗力の和がどうなるか？
     
 
 【Processingの説明】
    * Javaベースのプログラミング言語
    * 初心者でも実行結果をすぐに確認できる
    * アートやインタラクティブな表現を簡単に実装できる
    * 通常のプログラミング言語と異なり、main関数は不要
    * プログラムは次の2つの関数を中心に動作する：
      - setup(): プログラム起動時に一度だけ実行される初期化関数
      - draw(): フレームごとに繰り返し実行される関数（デフォルトで毎秒60回）
    * 描画用のウィンドウは size(幅のピクセル数, 高さのピクセル数) で設定する
    * この関数は setup() 内で一度だけ呼ぶのが基本
    * 自動的に width（幅のピクセル数）と height（高さのピクセル数）が設定される
      - widthとheightは予約語のため、同名の変数は宣言不可
    * 豊富な描画機能（線、図形、色指定、画像処理など）が標準で組み込まれている
    * インタラクション（マウス、キーボード入力）も簡単に扱える
    * 一般的な数学関数（三角関数、乱数生成など）も組み込み済み
      - ただし浮動小数点数のデフォルトは 単精度（float）
      - 倍精度（double）を使いたい場合は、JavaのMathクラスを使う
     
 【実行方法】 
    * このウィンドウの左上にある三角ボタンをクリックすると新しいウィンドウが表示される
    * そのウィンドウ上でマウスをクリックするとシミュレーションが開始する
    * もう一度クリックすると再び一時停止
 
 */



// -----------------
//    数学定数
// -----------------
final double TWOPI = 2 * Math.PI;
               // 円周率 π の2倍 


// -----------------
//    制御パラメータ
// -----------------
final double PARAM_TIME_STEP = 1.e-3;  
               // 時間刻み幅 dt の制御係数（1e-2は大きすぎ。1e-3程度以下。）
final int    PARAM_VIEW_SPEED = 2000;   
               // 表示速度の加速係数（大きいほど高速再生に見える）
final double PARAM_SPRING_DAMPER =1.0 ;      
//final double PARAM_SPRING_DAMPER =10.0;      
               // 砂粒のダンパーの減衰係数。1.0なら臨界減衰率               
final double PARAM_AVERAGE_TIME_SPAN_HINT = 1.0; 
               // 時間変動するデータの平均値をとる時間（単位は秒）の目安               
final double PARAM_HOURGLASS_TIME_IN_SECOND = 10.0;    //単位は秒(？)
//final double PARAM_HOURGLASS_TIME_IN_SECOND = 20;


// -----------------
//    基本整数
// -----------------
final int NSGIP = 100;      
            // Number of Sand Grains In a Pillar 一本の砂柱の中の砂粒の数
final int NSP = 3; 
            //Number of Sand Pillars　砂柱の数

// --------------
//    重力加速度
// --------------
final double GRAVITY_ACCELERATION = 9.80665;   
               // 重力加速度 (m/s^2)


// -----------
//    砂時計
// -----------
final double HOURGLASS_HEIGHT = 1.0;  
               // 砂時計の高さ (m)
final double HOURGLASS_HALF_HEIGHT = HOURGLASS_HEIGHT / 2;  
               // 半分
final double HOURGLASS_SAND_GRAIN_RELEASE_SECOND = PARAM_HOURGLASS_TIME_IN_SECOND / NSGIP;


// --------
//    砂
// --------
final double SAND_TOTAL_MASS = 0.1;  
               // 砂全体の質量（kg)
final double SAND_GRAIN_MASS = SAND_TOTAL_MASS / NSGIP;                     
               // 砂粒一つの質量 (kg)
final double SAND_PILE_HEIGHT = HOURGLASS_HALF_HEIGHT * 0.8;
               // 砂粒層の厚さ (m)
final double SAND_GRAIN_DIAMETER = SAND_PILE_HEIGHT / NSGIP;                    
               // 砂粒の直径 (m) 
final double SAND_GRAIN_RADIUS = SAND_GRAIN_DIAMETER / 2;  
               // 砂粒の半径 (m)


// -----------------------------------------------
//    自由落下速度（速度のスケールの評価に使う）
// -----------------------------------------------
final double FREE_FALL_VELOCITY 
               = Math.sqrt(2 * GRAVITY_ACCELERATION * HOURGLASS_HALF_HEIGHT);
                  // 自由落下速度の最大値 (m/s)
final double SMALL_TIME_SCALE_BY_VELOCITY 
               = SAND_GRAIN_RADIUS / FREE_FALL_VELOCITY;
                   // 自由落下速度が砂粒の半径を横切るのにかかる時間


// -----------------
//    計算対象の範囲
// -----------------
final double SIMULATION_REGION_Y_MAX =  HOURGLASS_HALF_HEIGHT;
final double SIMULATION_REGION_Y_MIN = -HOURGLASS_HALF_HEIGHT;
final double SIMULATION_REGION_X_MAX 
               = SIMULATION_REGION_Y_MAX * ((double)width / (double)height);
final double SIMULATION_REGION_X_MIN 
               = SIMULATION_REGION_Y_MIN * ((double)width / (double)height);


// --------------------
//    砂粒の弾性（バネ）
// --------------------
     /*
           落下する砂粒が床面に到達したら跳ね返るようにするため、砂粒の
           弾性としてバネ・ダンパーモデルを採用する。
           
           これは以下のようなモデルである。床面にぶつかった砂粒はわずかに
           変形するはずである。その歪みに比例した強さでもとの球形に戻る
           ような反発力が作用する。その反発力をバネで表現する。つまり
           砂粒の中心と床面との距離 d を常に測っておき、d が砂粒の
           半径 a よりも小さくなったら、その差 (a-d) に比例する力が
           作用して、床が砂粒を押し返すと仮定する。距離に比例する
           力はバネで表現できる。ただし、ただのバネを想定すると、砂粒が
           床面に何度衝突してもエネルギーの損失がないので、床に跳
           ね返された砂粒はもとの速度で上空に投げ返されてしまう。
           
           実際の砂時計の砂粒はもちろん床面に（あるいは先に落ちて床の上に
           積み重なった砂の層の上面に）衝突する過程でほとんど一瞬で
           運動エネルギーを失うので、上空に跳ね返されることはないで
           あろう。そこで、衝突の際のエネルギー損失をダンパーでモデル
           する。ダンパーとは砂粒の変形速度に比例し、速度の向きとは
           反対方向に作用する力である。速さに比例する空気抵抗のモデル
           と同じといえる。
           
           バネ・ダンパーモデルのバネ定数 k をどう決めるかは自由であるが、
           ここでは以下の方針で決める。
           
           kの値が大きいほど硬い砂粒、つまり実際の砂粒を模擬することになるが、
           あまり大きすぎるとバネ振動の時間スケールで決まる時間ステップが
           小さくなりシミュレーションの進行が遅くなる。そこでここでは、
           
             バネ振動の周期 = 砂粒の自由落下速度で決まる小さい時間スケール
             
           という条件で k の値を決めることにする。右辺の「砂粒の自由落下
           速度で決まる小さい時間スケール」というのは、砂粒が砂時計の中を
           落下するときの最大速度（つまり一番下の床面まで自由落下したとき
           の速度）が砂粒の半径を横切るのにかかる時間である。砂粒と床面
           との衝突を検知するためにはこの時間と比較して十分小さい時間間隔
           で砂粒が床面に到達したかどうか判断（衝突検知）する必要があるので、
           この「砂粒の自由落下速度で決まる小さい時間スケール」はシミュレー
           ションの時間間隔を決めるいちばん厳しい条件である。
                      
           バネ・ダンパーモデルにおいてダンパーの強さ（減衰率）も自由に
           決められるパラメータである。減衰率が大きいと振動はしないが砂粒は
           なかなか変形しないし、減衰率が小さいとエネルギーの損失が小さいために
           バネ振動が長時間続く。最も効率的に振動を減衰させる減衰率は
           臨界減衰率と呼ばれ、バネ定数 k と砂粒の質量 m だけで決まる。
           ここではこの臨界減衰率を基準としてその何倍か、
           と方法で砂粒のダンパーの減衰率を指定する。
             
      */
final double SPRING_CONST   // バネ定数 (N/m)
               = SAND_GRAIN_MASS * Math.pow(TWOPI / SMALL_TIME_SCALE_BY_VELOCITY, 2);
                    /*
                        砂粒の質量をmとし、バネ定数をkとすると、バネの角振動数は 
                           omega = sqrt(k/m)
                        で、その周期は 
                           tau = TWOPI / omega = TWOPI / sqrt(k/m) 
                        である。いま、バネ定数 k を 
                           tau = SMALL_TIME_SCALE_BY_VELOCITY
                        という条件から決める。つまり
                           TWOPI / sqrt(k/m) = SMALL_TIME_SCALE_BY_VELOCITY
                        である。両辺を2乗して変形すると 
                          k = m * (TWOPI / SMALL_TIME_SCALE_BY_VELOCITY)^2
                    */
final double SPRING_FREQUENCY_OMEGA = Math.sqrt(SPRING_CONST / SAND_GRAIN_MASS);
               // バネ振動の角振動数 (1/s)
final double SPRING_PERIOD = TWOPI / SPRING_FREQUENCY_OMEGA;
               // バネ振動の周期 (s)
final double DAMPER_CONST_CRITICAL_VALUE 
               = 2 * Math.sqrt(SAND_GRAIN_MASS * SPRING_CONST);
                  // バネ・ダンパーモデルにおける臨界減衰係数
final double DAMPER_CONST 
               = PARAM_SPRING_DAMPER * DAMPER_CONST_CRITICAL_VALUE;
                  // バネ・ダンパーモデルにおける減衰係数


// ------------------------------------------------
//    砂粒が床面と接触したかどうかの判定基準距離
// ------------------------------------------------
final double CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR = SAND_GRAIN_RADIUS;
           /*
            
                砂粒の中心と床面との距離を常に測っておき、その距離が
                   CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR
                よりも小さくなったら、砂粒が床面に接触していると判断する。
           
                Falling sand grain
                    
                       o o
                     o     o             o o
                    o   +   o          o     o     
                     o     o          o   +   o        --.     
                       o o             o     o             Distance between  
                   ___________       ____o_o____       __. sand grain and floorLower
                      Floor            Just touch 
                                         on the floorLower  
            */



// ----------------------------------------
//    砂柱毎の時間管理
// ----------------------------------------

Laptimes laptimes;

class Laptimes
{
  double[] laptime = new double[NSP];
  
  Laptimes()
  {
    randamize_laptime();
  }
  
  void randamize_laptime() {
    for (int p=0; p<NSP; p++) {
      laptime[p] = HOURGLASS_SAND_GRAIN_RELEASE_SECOND * Math.random(); //random()は0.0から1.0未満のdouble型の値を返す。
    }
  }

  double get_laptime(int p)
  {
    return laptime[p];
  }
  
  void reset_laptime(int p)
  {
    laptime[p] = sim.time;
    println("resetting laptime for p = " + p);
  }
  
}

// ----------------------------------------
//    シミュレーションの時刻とステップ数
// ----------------------------------------
Simulation sim;

class Simulation
{
  double time;
  int    nstep;
  double dt;
  
  boolean time_keeping_on;  // 計時しているかどうか
  
  String str_time;
  String str_nstep;
  String str_dt;
  
  Simulation() 
  {
    time = 0.0;  // シミュレーションの時刻 (s) 
    nstep = 0;   // シミュレーションのステップ数
    dt = PARAM_TIME_STEP * Math.min(SMALL_TIME_SCALE_BY_VELOCITY, 
                                    SPRING_PERIOD);
             // バネの振動と自由落下速度のそれぞれの時間スケールのうち小さい方を
             // 採用し、その定数倍を時間刻み幅 dt とする
             // PARAM_TIME_STEP は10^{-3} 程度がいい。
    str_time  = nfs((float)time,1,8);
    str_nstep = nfs(nstep, 9);
    str_dt    = nfs((float)dt,1,8);
    time_keeping_on = false;
  }
  
  void start_time_keeping()
  {
    time_keeping_on = true;
  }
  
  void stepIncrement()
  {
    nstep += 1;
    str_nstep = nfs(nstep, 9);
  }
  
  void timeIncrement()
  {
    time += dt;
    str_time  = nfs((float)time,1,8);
  }
}


// ----------------------------------------------------
//    シミュレーションの実行と一時停止用のトグル変数
// ----------------------------------------------------
boolean RunningStateToggle = false;  
          // false なら一時停止状態。マウスクリックで切り替え。


// ----------------
//    砂時計の床
// ----------------
Floor floorLower;

class Floor
{
  double level_y;            
           // 床面のy座標
  int[] touching_grain = new int[NSP];     
           // 一本の砂柱の中の、どの砂粒と接触する可能性があるか
  double normal_force;       
           // 接触している砂粒から床が受けるバネの力の反作用（=垂直抗力）
  double draw_width;         
           // 床面を長方形で表示するときの幅 (m)
  double draw_height;        
           // 床面を長方形で表示するときの高さ (m)
  double draw_width_left_x;  
           // 床面長方形左端のx座標 (m) 
  
  
 
    Floor(int touching_grain, double level_y)
    {
      for(int i=0; i < NSP; i++){
        this.touching_grain[i] = touching_grain;    // 床と相互作用する粒子番号
      }
      
      this.level_y = level_y;
      draw_width = ( SIMULATION_REGION_X_MAX 
                   - SIMULATION_REGION_X_MIN ) * 0.75;
      draw_height = SAND_GRAIN_RADIUS*4;
      draw_width_left_x = - draw_width / 2;
      
    }
  
  void resetNormalForce()
  {
    normal_force = 0.0;
  }
  
  double getNormalForce()
  {
    return normal_force;
  }
  
  void switch_touching_grain(int i)
  {
    if( this.touching_grain[i] < NSGIP ) {
      this.touching_grain[i] += 1;
    }
    
  }
  
  void draw() 
  {    
    stroke(0);
   
    fill(180);
    rect(mapx(draw_width_left_x), 
         mapy(level_y - draw_height),
         mapx(draw_width), 
         mapy(draw_height));
  }    
}


// -----------------------------
//   上の床面（オリフィスのある面
// -----------------------------

Floor floorUpper;

 
// --------------------
//    砂粒とその配列
// --------------------

Grain[][] grains = new Grain[NSP][NSGIP];     // [砂柱の数][砂粒の数]
double[] phaseShiftForEachPillar = new double[NSP];     //それぞれの砂柱について最初のドロップがあるまでの時間を格納。
  

  /*
          When NSGIP=3
      
              y   
                      
              |       sand(i=2)
              |    .   
              | . 
              o       sand(i=1)
              |     .   
        ------+---.------------- x
              | . 
              o       sand(i=0)
              |     .    
              |   .     Floor
              | .    .   
           ___o___.            
  
   */

class Grain
{
  double pos_x;  // 砂粒のx座標（時間変化しない）
  double pos_y;  // 砂粒のy座標（時間変化する）
  double vel_y;  // 砂粒の速度のy成分 （時間変化する）
  
  
  
  Grain(double x, double y, double vy)
  {
    pos_x = x;
    pos_y = y;
    vel_y = vy;
}
  
  
  
  void draw()
  {
    stroke(0);
    fill(255, 255, 190);
    circle(mapx(pos_x),                  // x座標 
           mapy(pos_y),                  // y座標
           mapy(SAND_GRAIN_RADIUS*2));   // 円の直径
  }  
} 


// ---------------------
//    データ解析：平均値
// ---------------------

class Average
{
  /*
  
    連続して得られる多数の数値の平均値を計算する。
    
    全数値の平均値ではなく、一番最近の N 個のデータの平均をとる。
    
    たとえば、 N=10^6(=1000000) のように大きい場合、平均値 
      v_mean = (v_0 + v_1 + ... + v_999999) / 1000000 
    を求めるのに、大きさ1000000の配列を用意して計算するのは
    メモリの無駄である。
    
    ここでは「平均の平均は全体の平均である」という事実を使って以下のように
    して v_mean を求める。
    
    簡単のため N=4 のときで説明すると
        v_mean = ( v_0 + v_1 + v_2 + v_3 ) / 4
    を求めるのに、
        v_mean01 = ( v_0 + v_1 ) / 2
        v_mean23 = ( v_2 + v_3 ) / 2
    をまずは計算して
        v_mean = ( v_mean01 + v_mean23 ) / 2
    とする。
    
    同様に、N=100 のときには      
    
    v_mean = ( 
                v_0 + v_1 + ... + v_99
              ) 
                / 100
           = (   
                 ( v_0   + v_1   + ... + v_9  ) / 10 
               + ( v_10 + v_11 + ... + v_19 ) / 10
               + ( v_20 + v_21 + ... + v_29 ) / 10
               .
               .
               + ( v_90 + v_91 + ... + v_99 ) / 10
             ) 
               / 10
    
    として v_mean を計算する。
              
   */

  int       size_for_work_array;  // 上のコメントで N に相当するもの
  double    time_span_for_work_array;
  double    actual_time_span_for_average;
  double[]  array_of_rawdata;
  int       ctr_for_array_of_rawdata;  // ctr = "counter"
  double    avrg_of_array_of_rawdata;
  String    str_avrg_of_array_of_rawdata;
  double[]  array_of_mean_of_array_of_rawdata;
  double[]  array_of_avrg_of_array_of_rawdata;
  int       ctr_for_array_of_avrg_of_array_of_rawdata;
  double    avrg_of_array_of_avrg_of_array_of_rawdata;
  String    str_avrg_of_array_of_avrg_of_array_of_rawdata;
  double[]  array_of_avrg_of_array_of_avrg;
  int       ctr_for_array_of_avrg_of_array_of_avrg;


  double take_average(double[] array)
  {
    double sum = 0.0;
    for (int i=0; i<size_for_work_array; i++) {
      sum += array[i];
    }
    return sum / size_for_work_array;
  }

  Average(double dt)
  {
    size_for_work_array = (int)Math.round(Math.sqrt(PARAM_AVERAGE_TIME_SPAN_HINT/dt));
    time_span_for_work_array = dt * size_for_work_array;
    actual_time_span_for_average = dt * Math.pow(size_for_work_array,2);

    println("<Average>          size_for_work_array = ", size_for_work_array );
    println("<Average> time_span_for_work_array     = ", time_span_for_work_array );
    println("<Average> actual_time_span_for_average = ", actual_time_span_for_average );

    array_of_rawdata                  = new double[size_for_work_array];
    array_of_avrg_of_array_of_rawdata = new double[size_for_work_array];
    array_of_avrg_of_array_of_avrg    = new double[size_for_work_array];
      
    ctr_for_array_of_rawdata                  = 0;
    ctr_for_array_of_avrg_of_array_of_rawdata = 0; 
    ctr_for_array_of_avrg_of_array_of_avrg    = 0; 
    
    str_avrg_of_array_of_rawdata                  = " -.--------";
    str_avrg_of_array_of_avrg_of_array_of_rawdata = " -.--------";
  
    for (int n=0; n<size_for_work_array; n++) {
      array_of_rawdata[n]                  = 0.0;      
      array_of_avrg_of_array_of_rawdata[n] = 0.0;
      array_of_avrg_of_array_of_avrg[n]    = 0.0;
    }
  }

  void register(double rawdata)
  {    
      
    /*  
       When size_for_work_array = 4
       
       1 3 2 0 3 2 4 5 5 4 7 6 8 7 6 9  <== array_of_rawdata
        \___/   \___/   \___/   \___/   
          3       7      11      15     <== array_of_avrg_of_array_of_rawdata
           \____________________/    
                      9                 <== avrg_of_array_of_avrg_of_array_of_rawdata

     */

    array_of_rawdata
        [ctr_for_array_of_rawdata] = rawdata;
         ctr_for_array_of_rawdata += 1;    
    if ( ctr_for_array_of_rawdata == size_for_work_array ) {
         ctr_for_array_of_rawdata = 0; // reset      
          avrg_of_array_of_rawdata = take_average(array_of_rawdata);
      str_avrg_of_array_of_rawdata = nfs((float)avrg_of_array_of_rawdata,1,8);
      array_of_avrg_of_array_of_rawdata
          [ctr_for_array_of_avrg_of_array_of_rawdata] = avrg_of_array_of_rawdata;
           ctr_for_array_of_avrg_of_array_of_rawdata += 1;      
      if ( ctr_for_array_of_avrg_of_array_of_rawdata == size_for_work_array ) {
           ctr_for_array_of_avrg_of_array_of_rawdata = 0; // reset
          avrg_of_array_of_avrg_of_array_of_rawdata = take_average(array_of_avrg_of_array_of_rawdata);
      str_avrg_of_array_of_avrg_of_array_of_rawdata = nfs((float)avrg_of_array_of_avrg_of_array_of_rawdata,1,8);
        array_of_avrg_of_array_of_avrg
            [ctr_for_array_of_avrg_of_array_of_avrg] = avrg_of_array_of_avrg_of_array_of_rawdata;
             ctr_for_array_of_avrg_of_array_of_avrg += 1;      
        if ( ctr_for_array_of_avrg_of_array_of_avrg == size_for_work_array ) {
          ctr_for_array_of_avrg_of_array_of_avrg = 0; // reset
        }
      }
    }
  }
  
  void drawSimplePlot()
  {
    double plot_x_min = 0.35*width;
    double plot_x_max = 0.95*width;
    double plot_y_min = 0.85*height;
    double plot_y_max = 0.95*height;
    double plot_width  = plot_x_max - plot_x_min;
    double plot_height = plot_y_max - plot_y_min;
    double weight_min = SAND_GRAIN_MASS*NSGIP * ( 1 - 5.e-1 );
    double weight_max = SAND_GRAIN_MASS*NSGIP * ( 1 + 5.e-1 );
    double plot_dx = plot_width / size_for_work_array;
    double plot_dy = plot_height / ( weight_max - weight_min );

    double y_minus = plot_y_min;
    double y_zero  = plot_y_min + plot_dy*(SAND_GRAIN_MASS*NSGIP-weight_min);
    double y_plus  = plot_y_min + plot_dy*(weight_max            -weight_min);

    float fx_min   = (float)plot_x_min;
    float fx_max   = (float)plot_x_max;
    float fy_minus = (float)y_minus;
    float fy_zero  = (float)y_zero;
    float fy_plus  = (float)y_plus;
    
    translate(0, height);
    scale(1, -1); 

    strokeWeight(0.5);
    stroke(0,255,0);
    line( fx_min, fy_zero, fx_max, fy_zero);
    stroke(255,0,0);
    line( fx_min, fy_plus, fx_max, fy_plus);
    stroke(0,0,255);
    line( fx_min, fy_minus, fx_max, fy_minus);

    stroke(0);

    fill(100);
    for (int i=0; i<size_for_work_array; i++) {
      double plot_x = plot_x_min + plot_dx * i;
      double weight = array_of_avrg_of_array_of_rawdata[i];
      if ( weight >= weight_min && weight <= weight_max ) {
        double plot_y = plot_y_min + plot_dy * (weight - weight_min);
        circle((float)plot_x, (float)plot_y, 2);
      }
    }
    
    fill(255,0,0);
    for (int i=0; i<size_for_work_array; i++) {
      double plot_x = plot_x_min + plot_dx * i;
      double weight = array_of_avrg_of_array_of_avrg[i];
      if ( weight >= weight_min && weight <= weight_max ) {
        double plot_y = plot_y_min + plot_dy * (weight - weight_min);
        circle((float)plot_x, (float)plot_y, 4);
      }
    }
  }

} 



// ----------------------------------
//    データ解析：エネルギーの計算
// ----------------------------------

class Energy
{
  String str_total_energy;
  
  Energy()
  {
    str_total_energy = " -.--------";
  }
  
  double getKineticEnergy()
  {
    /*
          砂粒の運動エネルギーの和   \sum_{i=0}^{NSGIP-1} (1/2) * m * v_i^2
     */
    double sum = 0.0;
    
    for(int p=0; p < NSP; p++){
      for(int i=0; i<NSGIP;i++){
        double vel_y = grains[p][i].vel_y;
        double vel_y_sq = vel_y * vel_y;
        sum += 0.5*SAND_GRAIN_MASS*vel_y_sq;
      }
    }
    return sum;
  }
  
  double getGravityPotential()
  {
    /*
     
         砂粒の重力ポテンシャルの和   \sum_{i=0}^{NSGIP-1} m * g * h_i
         
         高さh_i の基準は砂時計の一番下のy座標にする。
         
         +-----+
         |     |
         |     |
         |     |
          \   /
           \ /
           / \
          / o-\- - - - - +  
         |  .  |         |
         |  .  |    height for gravity potential
         |  .  |         |
         +--+--+ - - - - +
         
     */    

    double sum = 0.0;
    
    
    for (int p=0; p<NSP; p++){
      for(int i=0; i<NSGIP; i++){
        double height_of_grain = grains[p][i].pos_y - floorLower.level_y;
        sum +=  SAND_GRAIN_MASS 
            * GRAVITY_ACCELERATION 
            * height_of_grain;
      }
    }
    return sum;
  }
  
  double getSpringPotential()
  {
    /*
          弾性ポテンシャルの和 U = U_s + U_f
          
          弾性ポテンシャルは2種類ある。

              (1) 砂粒と砂粒の間の相互作用のポテンシャル

                     今の場合、仮想バネを想定しているのでバネのポテンシャル
                        U_s = (1/2) * k' * L^2
                     ここで k' は仮想バネのバネ定数で、
                     Lは隣り合う二つの質点の間に想定している仮想バネの伸び、
                     つまり二つの質点の間の距離と仮想バネの自然長との差である。

              (2) 砂粒と床面が接触しているときのポテンシャルエネルギー
                     
                     砂粒の中心と床面との距離を d とすると、dがあらかじめ
                     設定した長さ
                       CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR
                     よりも短いときに床面との接触による力が生じる。
                     逆に d > CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR
                     のときには接触していないので、力は生じない。
                     
                     つまり                     
                      L = CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR - d
                     が正のときにだけ力が働く。
                     
                     砂粒の弾性バネ定数を k, 自然長を L0 とすると
                     接触時に作用する力の大きさは
                       f = k * (L-L0)
                     である。ポテンシャルは
                       U_f = (1/2) * k * (L-L0)^2 
                     である。
     
     */
    double sum = 0.0;

    for(int p=0; p<NSP;p++){
      for(int i=0; i<NSGIP; i++){
        
        if (i < NSGIP-1){  //上の粒子と繋がってるバネのエネルギー。j = NSGIP-1 は定義されない
          double dist_to_upper_neighbor = grains[p][i+1].pos_y - grains[p][i].pos_y;
          positive_check( dist_to_upper_neighbor, "dist_to_upper_neighbor < 0?" );
      
          double overlap_upper = SAND_GRAIN_DIAMETER - dist_to_upper_neighbor;   //砂粒同士の接触判定は距離が直径以下どうかで判断
      
          if(overlap_upper > 0){
            double overlap_upper_sq = Math.pow(overlap_upper,2);
            sum += 0.5*SPRING_CONST*overlap_upper_sq;
          }
        }
        
        if (i > 0 ){    //下の粒子と繋がってるバネのエネルギー。j = 0 は定義されない
          double dist_to_lower_neighbor = grains[p][i].pos_y - grains[p][i-1].pos_y;
          positive_check( dist_to_lower_neighbor, "dist_to_upper_neighbor < 0?" );
      
          double overlap_lower = SAND_GRAIN_DIAMETER - dist_to_lower_neighbor;   //砂粒同士の接触判定は距離が直径以下どうかで判断
      
          if(overlap_lower > 0){
            double overlap_lower_sq = Math.pow(overlap_lower,2);
            sum += 0.5*SPRING_CONST*overlap_lower_sq;
          }
        }
        
        if ( i==floorLower.touching_grain[p] ) { 
          double dist_to_floorLower = grains[p][i].pos_y - floorLower.level_y;
          positive_check( dist_to_floorLower, "dist_to_floorLower < 0?" );
          double overlap_floorLower = CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR 
                         - dist_to_floorLower;
          if ( overlap_floorLower> 0 ) {
            double overlap_floorLower_sq = overlap_floorLower * overlap_floorLower;
            sum += 0.5*SPRING_CONST*overlap_floorLower_sq;
          }
        }
        
        if ( i==floorUpper.touching_grain[p] ) { 
          double dist_to_floorUpper = grains[p][i].pos_y - floorUpper.level_y;
          positive_check( dist_to_floorUpper, "dist_to_floorLower < 0?" );
          double overlap_floorUpper = CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR 
                         - dist_to_floorUpper;
          if ( overlap_floorUpper> 0 ) {
            double overlap_floorUpper_sq =overlap_floorUpper * overlap_floorUpper;
            sum += 0.5*SPRING_CONST*overlap_floorUpper_sq;
          }
        }
        
      }
    }
    
    
    return sum;
  }
  
  double getTotalEnergy()
  {
    double kinetic   = getKineticEnergy();
    double gravity   = getGravityPotential();
    double potential = getSpringPotential();
    double total_energy = potential + kinetic + gravity ;
    
    str_total_energy = nfs((float)total_energy,1,8);
    return total_energy;
  }
} 


// --------------------------
//    データ解析機能の集成
// --------------------------

Analyser analyser;

class Analyser
{
  Average average;
  Energy  energy;
  Analyser(double dt) 
  {
    average = new Average(dt);
    energy  = new Energy();
  }
}


// ----------------------------------------
//    最初に一度だけ自動的に呼ばれる関数
// ----------------------------------------

void setup() 
{
  size(500, 800);    // 実行画面を表示するウィンドウの横幅(width)と縦幅（height）の設定
  background(255);   // 背景色の指定
  initialize();      // シミュレーションの初期設定
  
  println("PARAM_TIME_STEP = ", PARAM_TIME_STEP);
  println("SMALL_TIME_SCALE_BY_VELOCITY = ", SMALL_TIME_SCALE_BY_VELOCITY);
  println(" SPRING_PERIOD   = ", SPRING_PERIOD);
  println("dt   = ", sim.dt);
  println( "     NSGIP = ", NSGIP);
}


// ----------------------------------
//    シミュレーションの初期化処理
// ----------------------------------
void initialize()
{
  /*
      setup関数の中から呼び出される
   */
  int touching_grain = 0;
  
  sim = new Simulation();
  floorLower = new Floor( touching_grain, SIMULATION_REGION_Y_MIN*0.9 );
  analyser = new Analyser(sim.dt);
  floorUpper = new Floor(touching_grain, - SIMULATION_REGION_Y_MIN*0.0);   //Upperfloorの設定
  laptimes = new Laptimes();
  
  double separation = SAND_GRAIN_DIAMETER;    // 砂粒の初期距離は直径でok
  
  
  for (int p=0; p<NSP; p++){
    double x = floorLower.draw_width_left_x +  floorLower.draw_width * p / (NSP-1) ; //植木算に注意　NSP-1でok
    for (int i=0; i<NSGIP; i++){
      double  y = - SIMULATION_REGION_Y_MIN*0.25 + separation*i;
      double vy = 0.0;
      grains[p][i] =new Grain(x,y,vy);
    }
  }
}


// -----------------------------------------------
//    アサート（実行時の問題発生を早期発見する）
// -----------------------------------------------
void positive_check(double must_be_positive, String last_will)
{
  /*
       常に正の値であるはずの第一引数 must_be_positive の値が
       本当に正になっているかどうかの確認。万が一、この値が正に
       なっていなかったら、第二引数で渡された文字列を出力してから
       プログラムを終了する。
   */
  if ( must_be_positive <= 0.0 ) {
    println("#Fatal error: " + last_will);
    exit();
  }
}


// -----------------------------------
//    運動方程式 equation of motion                            //どの砂柱について計算しているのかを知る必要がある
// -----------------------------------
void equationOfMotion(double  posy[],
                      double  vely[],                      
                      double dposy[],
                      double dvely[],
                      double sim_dt, 
                      int  pillar_id)   //砂柱の番号 (0 以上 NSP-1以下)
{
  /*
       一つの砂柱の中の全砂粒の運動をニュートンの運動方程式にしたがって解く
       rungeKutta4 関数の中から呼び出される
   */
  
  for (int i=0; i<NSGIP; i++) {   
    
    //-------------------------
    //  上の砂粒との相互作用 
    //-------------------------
    double spring_force_from_upper_neighbor = 0.0;  // default
    double damper_force_from_upper_neighbor = 0.0;
    
    if ( i < NSGIP-1 ) {
      double dist_to_upper_neighbor = posy[i+1] - posy[i];
      positive_check( dist_to_upper_neighbor, "dist_to_upper_neighbor < 0?" );
      
      double overlap_upper = SAND_GRAIN_DIAMETER - dist_to_upper_neighbor;    //砂粒同士の接触判定は距離が直径以下どうかで判断
      
      if (overlap_upper > 0){
      
      spring_force_from_upper_neighbor 
        =  SPRING_CONST * ( dist_to_upper_neighbor
                                      - SAND_GRAIN_DIAMETER);
      damper_force_from_upper_neighbor = - DAMPER_CONST * (vely[i]-vely[i+1]);
      }

        
    }

    //-------------------------
    //  下の砂粒との相互作用
    //-------------------------
    double spring_force_from_lower_neighbor = 0.0;
    double damper_force_from_lower_neighbor = 0.0;

    if ( i > 0 ) {
      double dist_to_lower_neighbor = posy[i] - posy[i-1];
      positive_check( dist_to_lower_neighbor, "dist_to_lower_neighbor < 0?" );
      
      double overlap_lower = SAND_GRAIN_DIAMETER - dist_to_lower_neighbor;
      
      if(overlap_lower > 0){
        spring_force_from_lower_neighbor 
        = - SPRING_CONST * ( dist_to_lower_neighbor
                                      - SAND_GRAIN_DIAMETER);
        damper_force_from_lower_neighbor = - DAMPER_CONST * (vely[i]-vely[i-1]);
      }
    
    }

    //---------------------------
    //  砂時計の底面からの抗力 
    //---------------------------
    double spring_force_from_floorLower = 0.0;
    double damper_force_from_floorLower = 0.0;    
    
    if ( i==floorLower.touching_grain[pillar_id] ) {
      double dist_to_floorLower = posy[i] - floorLower.level_y;
      positive_check( dist_to_floorLower, "dist_to_floorLowerLower < 0?" );
      
      double overlap = CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR 
                     - dist_to_floorLower;
      if ( overlap > 0 ) {
        spring_force_from_floorLower = SPRING_CONST*overlap; // 上向きの力なので正
        damper_force_from_floorLower = - DAMPER_CONST * vely[i];        
        floorLower.normal_force = spring_force_from_floorLower 
                           + damper_force_from_floorLower;
      }
    }
    
    
    //--------
    //  オリフィスからの抗力 
    //--------
    
    double spring_force_from_floorUpper = 0.0;
    double damper_force_from_floorUpper = 0.0;    
    
    if ( i==floorUpper.touching_grain[pillar_id] ) {
      double dist_to_floorUpper = posy[i] - floorUpper.level_y;
      positive_check( dist_to_floorUpper, "dist_to_floorLowerLower < 0?" );
      
      double overlap = CONTACT_DISTANCE_BETWEEN_GRAIN_AND_FLOOR 
                     - dist_to_floorUpper;
      if ( overlap > 0 ) {
        spring_force_from_floorUpper = SPRING_CONST*overlap; // 上向きの力なので正
        damper_force_from_floorUpper = - DAMPER_CONST * vely[i];        
        floorUpper.normal_force = spring_force_from_floorUpper 
                           + damper_force_from_floorUpper;
      }
    }
    

    //--------
    //  重力 
    //--------
    double gravity_force = - SAND_GRAIN_MASS*GRAVITY_ACCELERATION;
    

    //-----------------------
    //  全ての力の和をとる 
    //-----------------------
    double force_total = spring_force_from_lower_neighbor + damper_force_from_lower_neighbor   //damperの力を足してる
                       + spring_force_from_upper_neighbor + damper_force_from_upper_neighbor
                       + spring_force_from_floorLower
                       + damper_force_from_floorLower
                       + spring_force_from_floorUpper
                       + damper_force_from_floorUpper
                       + gravity_force;

    dposy[i] = vely[i] * sim_dt;                       // dy = vy * dt
    dvely[i] = force_total * sim_dt / SAND_GRAIN_MASS; // dvy = (fy/m)*dt
  }
}


// -----------------------------------------------------
//    メートル単位のx座標をピクセル単位の横座標に変換
// -----------------------------------------------------

float mapx(double x) 
{
  /*
                 (x,y) = physical unit coords. 
       (map(x),map(y)) = pixel coords.
   */
  double xmax = SIMULATION_REGION_X_MAX;
  double xmin = SIMULATION_REGION_X_MIN;
  double scale = width/(xmax-xmin);
  return map((float)x, (float)xmin, (float)xmax, 
             (float)(scale*xmin), (float)(scale*xmax));
}


// -----------------------------------------------------
//    メートル単位のy座標をピクセル単位の縦座標に変換
// -----------------------------------------------------
float mapy(double y) 
{
  /*
                 (x,y) = physical unit coords. 
       (map(x),map(y)) = pixel coords.
   */
  double ymax = SIMULATION_REGION_Y_MAX;
  double ymin = SIMULATION_REGION_Y_MIN;
  double scale = height/(ymax-ymin);
  return map((float)y, (float)ymin, (float)ymax, 
             (float)(scale*ymin), (float)(scale*ymax));
}



// ---------------------------------------------
//    Processingが1秒間に多数回自動的に呼び出す関数
// ---------------------------------------------
void draw() 
{
  /*
       デフォルトでは1秒間に60回
       シミュレーションの時間発展を進める関数 rungeKutta4() はここで呼ぶ       
   */
  background(255);
  
  push();  
    draw_sand_grains_and_floorLowers();
    draw_text_in_window();
  pop();

  push();
    analyser.average.drawSimplePlot();
  pop();

}

// -------------------
//    砂粒と床を表示
// -------------------
void draw_sand_grains_and_floorLowers()
{
  stroke(0, 0, 255);

  float ymax = (float)SIMULATION_REGION_Y_MAX;
  float ymin = (float)SIMULATION_REGION_Y_MIN;
  
  translate(width/2, height*ymax/(ymax-ymin));
  scale(1, -1);

  draw_grains();
  floorLower.draw();
  floorUpper.draw();                                  //オリフィスの描画

  
  if ( RunningStateToggle ) {
    for (int n=0; n<PARAM_VIEW_SPEED; n++) { // to speed up the display

      floorLower.resetNormalForce();                 // reset
      floorUpper.resetNormalForce();
      
      if ( sim.time_keeping_on ) {
        for (int p=0; p<NSP; p++) {
          if( sim.time - laptimes.get_laptime(p) > HOURGLASS_SAND_GRAIN_RELEASE_SECOND ) {
println(" switching p = " + p + " time = " + sim.time + "laptime = " + laptimes.get_laptime(p) 
          + " HOURGLASS_SAND_GRAIN_RELEASE_SECOND = " + HOURGLASS_SAND_GRAIN_RELEASE_SECOND);     
              floorUpper.switch_touching_grain(p);                       
              laptimes.reset_laptime(p);
            }
        }
      }

      rungeKutta4();
      double hourglass_weight = (floorLower.getNormalForce() + floorUpper.getNormalForce())/ GRAVITY_ACCELERATION; //オリフィスにかかる力を加えている
      
      analyser.average.register(hourglass_weight);
      analyser.energy.getTotalEnergy();
      
      if ( sim.nstep%10000 == 0 ) {
        
        println(  "#nstep", sim.str_nstep, 
                       "t", sim.str_time, 
                  "energy", analyser.energy.str_total_energy,
                     "raw", nfs((float)hourglass_weight,1,8));
        println("average1", analyser.average.str_avrg_of_array_of_rawdata,
                "average2", analyser.average.str_avrg_of_array_of_avrg_of_array_of_rawdata
               );
      }
    }
  }
}


// ---------------
//    砂粒を描く
// ---------------
void draw_grains() 
{
  stroke(50, 100, 200);

  for(int p=0; p<NSP; p++){
    for(int i=0; i<NSGIP; i++){
      grains[p][i].draw();
    }
  }
}


// --------------------------------
//    ウィンドウ上に文字列を表示
// --------------------------------
void draw_text_in_window() 
{
  /* 
       文字列を複数行描く
   */
  fill(0, 0, 0);
  scale(1, -1);
  textSize(12);

  
  float x = -width * 0.45;
  float y = -height * 0.45;
  float separation = 16;
        
  text("t=" + sim.str_time,                                  x, y);
  y += separation;
  text("nstep=" + sim.str_nstep,                             x, y);
  y += separation;
  text("energy=" + analyser.energy.str_total_energy,         x, y);
  y += separation;
  text("average1=" + analyser.average.str_avrg_of_array_of_rawdata, x, y);
  y += separation;
  text("average2=" + analyser.average.str_avrg_of_array_of_avrg_of_array_of_rawdata, x, y);
  
}



// ------------------------------------------------------
//    マウスボタンが押されたときに自動的に呼ばれる関数
// ------------------------------------------------------
void mousePressed() {
   RunningStateToggle = !RunningStateToggle;
}


void keyPressed() {
  if ( key=='s' ) {
    sim.start_time_keeping();
    
    laptimes.randamize_laptime();
    
    println("Time keeping.");
  }
}  
