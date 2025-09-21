
void rungeKutta4Advance(int num, 
                        double[] posy,
                        double[] vely,
                        double[] posy1,
                        double[] vely1,
                        double[] dposy,
                        double[] dvely,
                        double factor)
{
  for (int i=0; i<num; i++) {
    posy[i] = posy1[i] + factor*dposy[i];
    vely[i] = vely1[i] + factor*dvely[i];
  }
}


void rungeKutta4()
{
  final double ONE_SIXTH = 1.0/6.0;
  final double ONE_THIRD = 1.0/3.0;

  double[] posyprev = new double[NSGIP];
  double[] posywork = new double[NSGIP];
  double[]   dposy1 = new double[NSGIP];
  double[]   dposy2 = new double[NSGIP];
  double[]   dposy3 = new double[NSGIP];
  double[]   dposy4 = new double[NSGIP];
  double[] velyprev = new double[NSGIP];
  double[] velywork = new double[NSGIP];
  double[]   dvely1 = new double[NSGIP];
  double[]   dvely2 = new double[NSGIP];
  double[]   dvely3 = new double[NSGIP];
  double[]   dvely4 = new double[NSGIP];
  
  for (int i=0; i<NSGIP; i++) {
    posyprev[i] = grains[i].pos_y;
    velyprev[i] = grains[i].vel_y;
  }

  //step 1 
  equationOfMotion(posyprev,
                   velyprev,
                     dposy1,
                     dvely1,
                     sim.dt);
  rungeKutta4Advance(NSGIP,
                     posywork,
                     velywork,
                     posyprev,
                     velyprev,
                       dposy1,
                       dvely1,
                          0.5);                        
  
  //step 2
  equationOfMotion(posywork,
                   velywork,
                     dposy2,
                     dvely2,
                     sim.dt);
  rungeKutta4Advance(NSGIP,
                     posywork,
                     velywork,
                     posyprev,
                     velyprev,
                       dposy2,
                       dvely2,
                          0.5);
                          
  //step 3
  equationOfMotion(posywork,
                   velywork,
                     dposy3,
                     dvely3,
                     sim.dt);
  rungeKutta4Advance(NSGIP,
                     posywork,
                     velywork,
                     posyprev,
                     velyprev,
                       dposy3,
                       dvely3,
                          1.0);

  //step 4
  equationOfMotion(posywork,
                   velywork,
                     dposy4,
                     dvely4,
                     sim.dt);
  
  //the result
  for (int i=0; i<NSGIP; i++) { 
    posywork[i] = posyprev[i] + (
                           ONE_SIXTH*dposy1[i]
                         + ONE_THIRD*dposy2[i]
                         + ONE_THIRD*dposy3[i]
                         + ONE_SIXTH*dposy4[i] 
                         );
    velywork[i] = velyprev[i] + (
                           ONE_SIXTH*dvely1[i]
                         + ONE_THIRD*dvely2[i]
                         + ONE_THIRD*dvely3[i]
                         + ONE_SIXTH*dvely4[i] 
                         );
  }
  
  for (int i=0; i<NSGIP; i++) {
    grains[i].pos_y = posywork[i];
    grains[i].vel_y = velywork[i];
  }
  
  sim.timeIncrement();
  sim.stepIncrement(); 
}  
