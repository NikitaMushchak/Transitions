#include <time.h>
#include <boost/random.hpp>
#include "MD.h"
#include "phys_constants.h"

void MD::setTemperatureNormal() //устанавливается температура
{
  boost::random::mt19937 generator(time(0));
  boost::random::normal_distribution<> distribution(0.0,sqrt(IC.Ek*P_1d_M));
  boost::random::variate_generator<boost::random::mt19937&, boost::random::normal_distribution<> > distr(generator, distribution);
  int_fast32_t i;
  double koeff;
  for(i=0; i<N; ++i)
    {
        V[i].a[0] = distr();
        V[i].a[1] = distr();
        V[i].a[2] = distr();
        Vav += V[i];
        //std::cerr<<"VV "<<i<<" "<<ParticleVelocity[i].x<<" "<<ParticleVelocity[i].y<<"\n";
    }
    Vav *= _1d_N;
    //std::cerr<<"T "<<Vav.a[0]<<" "<<Vav.a[1]<<" "<<Vav.a[2]<<"\n";
    Ek = 0;
    for(i=0; i<N; ++i)
    {
        V[i] -= Vav;
        Ek += boost::qvm::mag_sqr(V[i]);
    }
    //std::cerr<<"T "<<Ek<<" "<<Ek0<<" "<<IC.Ek<<"\n";
    Ek *= 0.5*P_M*_1d_N;
    koeff = sqrt(IC.Ek/Ek);

    Ek = 0;
    for(i=0; i<N; ++i)
    {
        V[i] *= koeff;
        Ek += boost::qvm::mag_sqr(V[i]);
    }
    Ek *= 0.5*P_M;
    T = Ek*_1d_N*_1d_k_PhysConst;
    //std::cerr<<"Temperature "<<Ek<<" "<<T<<" "<<0.5*P_M*_1d_N*Ek<<"\n";
}
