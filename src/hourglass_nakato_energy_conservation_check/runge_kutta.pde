
void rungeKutta4Advance(double[] posvel,
                        double[] posvel1,
                        double[] dposvel,
                        double factor)
{
  for (int n=0; n<NPOSVEL; n++) {
    posvel[n] = posvel1[n] + factor*dposvel[n];
  }
}


void rungeKutta4()
{
  final double ONE_SIXTH = 1.0/6.0;
  final double ONE_THIRD = 1.0/3.0;

  double[] posvelprev = new double[NPOSVEL];
  double[] posvelwork = new double[NPOSVEL];
  double[]   dposvel1 = new double[NPOSVEL];
  double[]   dposvel2 = new double[NPOSVEL];
  double[]   dposvel3 = new double[NPOSVEL];
  double[]   dposvel4 = new double[NPOSVEL];
  
  posvelprev[0] = grain.pos_y;
  posvelprev[1] = grain.vel_y;

  //step 1 
  equationOfMotion(posvelprev,
                     dposvel1,
                     sim.dt);
  rungeKutta4Advance(posvelwork,
                     posvelprev,
                       dposvel1,
                            0.5);                        
  
  //step 2
  equationOfMotion(posvelwork,
                     dposvel2,
                     sim.dt);
  rungeKutta4Advance(posvelwork,
                     posvelprev,
                       dposvel2,
                            0.5);
                          
  //step 3
  equationOfMotion(posvelwork,
                     dposvel3,
                     sim.dt);
  rungeKutta4Advance(posvelwork,
                     posvelprev,
                       dposvel3,
                            1.0);

  //step 4
  equationOfMotion(posvelwork,
                     dposvel4,
                     sim.dt);
  
  //the result
  for (int n=0; n<2; n++) { 
    posvelwork[n] = posvelprev[n] + (
                           ONE_SIXTH*dposvel1[n]
                         + ONE_THIRD*dposvel2[n]
                         + ONE_THIRD*dposvel3[n]
                         + ONE_SIXTH*dposvel4[n] 
                         );
  }
  
  grain.pos_y = posvelwork[0];
  grain.vel_y = posvelwork[1];
  
  sim.timeIncrement();
  sim.stepIncrement(); 
}  
