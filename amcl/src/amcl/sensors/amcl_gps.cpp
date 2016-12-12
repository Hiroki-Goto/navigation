
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include "amcl_gps.h"
using namespace amcl;

AmclGPSSensor::AmclGPSSensor(){
}


bool AmclGPSSensor::gpsSensorUpdata(pf_t *pf, GPSSensorData *gps_data){
    pf_gps_update_sensor(pf, (pf_gps_sensor_model) GPSModel, gps_data);
    return true;
}

double AmclGPSSensor::normal_distribution(double diff, double dispersion){
    double result;

    result = 1/(dispersion*sqrt(2*M_PI)) * exp(-1*diff*diff /(2*dispersion*dispersion) );

    return result;
}

double  AmclGPSSensor::GPSModel(AmclGPSSensor *gps_data, pf_sample_set_t* set){
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
      dis_diff = sqrt( (pose.v[0]-gps_data->x)*(pose.v[0]-gps_data->x) +(pose.v[1]-gps_data->y)*(pose.v[1]-gps_data->y) );
      if(gps_data->position_covariance_type == 3){
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
