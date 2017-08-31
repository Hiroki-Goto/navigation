
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include "amcl_gnss.h"
#include "pf_pdf.h"
using namespace amcl;

AmclgnssSensor::AmclgnssSensor(){
}


bool AmclgnssSensor::gnssSensorUpdata(pf_t *pf, gnssSensorData *gnss_data){
    pf_gnss_update_sensor(pf, (pf_gnss_sensor_model) gnssModel, gnss_data);
    return true;
}

double AmclgnssSensor::normal_distribution(double diff, double dispersion){
    double result;

    result = 1/(dispersion*sqrt(2*M_PI)) * exp(-1*diff*diff /(2*dispersion*dispersion) );

    return result;
}

double  AmclgnssSensor::gnssModel(AmclgnssSensor *gnss_data, pf_sample_set_t* set){
    double p;
    double total_weight;
    double dis_diff;
    pf_sample_t *sample;
    pf_vector_t pose;

    double float_disp=1.2;
    total_weight = 0.0;

    // Compute the sample weights
    for (int i = 0; i < set->sample_count; i++){
      sample = set->samples + i;
      pose = sample->pose;

      p = 1.0;
      dis_diff = sqrt( (pose.v[0]-gnss_data->x)*(pose.v[0]-gnss_data->x) +(pose.v[1]-gnss_data->y)*(pose.v[1]-gnss_data->y) );
      if(gnss_data->position_covariance_type == 3){
          //p += normal_distribution(dis_diff,0.2);
          p *= (1/(0.3*sqrt(2*M_PI)));
      }else{
          if(dis_diff <= float_disp){
              p *= (1/(float_disp*sqrt(2*M_PI)));
          }else{
              p *= normal_distribution(dis_diff,float_disp);
          }

      }
      sample->weight *= p;
      total_weight += sample->weight;
    }
    return (total_weight);
}

bool AmclgnssSensor::gnssSensor_reseting(pf_t *pf, gnssSensorData *gnss_data, double sample_num){

    pf_sample_set_t *set;
    pf_sample_t *sample;
    int particle_num;
    pf_vector_t pose;
    const double gnss_sigma = 2.0;

    set = pf->sets + pf->current_set;
    //ばら撒くパーティクルの個数
    int gnss_sample = set->sample_count * sample_num;
    std::cout << gnss_sample << std::endl;
     for(int i=0; i< gnss_sample; i++){
         particle_num = rand() % set->sample_count + 1;
         sample = set->samples + particle_num;
         pose = sample->pose;

         sample->pose.v[0] = gnss_data->x + pf_ran_gaussian(gnss_sigma);
         sample->pose.v[1] = gnss_data->y + pf_ran_gaussian(gnss_sigma);
         sample->pose.v[2] = pf_ran_gaussian(3.14);

     }
   return true;
}
