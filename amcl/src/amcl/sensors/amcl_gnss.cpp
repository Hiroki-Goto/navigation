
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

//行列演算用ライブラリ
#include <Eigen/Core>
#include <Eigen/LU> //逆行列と行列式がcoreにないためインクルード

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
    //std::cout << gnss_sample << std::endl;
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

double AmclgnssSensor::gnssPfKLD(pf_vector_t pf, gnssSensorData *gnss_data){
    double d = 2.0;   //対角成分の数

    double pf_sigma = 1;
    Eigen::Vector2d pf_position;
    Eigen::Matrix2d pf_sigma_mx;
    //std::cout << pf.v[0] << "\t" << pf.v[1] << std::endl;

    double gnss_sigma = 5.0;
    Eigen::Vector2d gnss_position;
    Eigen::Matrix2d gnss_sigma_mx;

    pf_position << pf.v[0], pf.v[1];
    pf_sigma_mx << pf_sigma, 0,
                   0, pf_sigma;

    //std::cout << gnss_data->x << "\t" << gnss_data->y << std::endl;
    gnss_position << gnss_data->x, gnss_data->y;
    gnss_sigma_mx << gnss_sigma, 0,
                  0, gnss_sigma;

    double kld = ( log(gnss_sigma_mx.determinant() / pf_sigma_mx.determinant())
                    + (gnss_sigma_mx.inverse()*pf_sigma_mx).trace()
                        + (pf_position-gnss_position).transpose()*gnss_sigma_mx.inverse()*(pf_position-gnss_position) - d) /2;
    //kld = log10(gnss_sigma_mx.determinant() / pf_sigma_mx.determinant());
    //std::cout << gnss_sigma_mx(0,0) << "\t" <<  gnss_sigma_mx(1,1) << std::endl;
    //std::cout << kld << std::endl;
    return kld;
}

void AmclgnssSensor::er(pf_t *pf, analysis_t *data){
    int reset_count = 3;
    pf_sample_set_t *set;
    pf_sample_t *sample;
    set = pf->sets + pf->current_set;
    //std::cout << data->dispersion.v[0] << data->dispersion.v[1] << std::endl;
    int reset_limit = ((int)data->dispersion.v[0] + (int)data->dispersion.v[1]) /2;
    //std::cout << reset_limit << std::endl;
    if(reset_count >= reset_limit){
        for(int i=0; i<set->sample_count; i++ ){
            sample = set->samples + i;
            ///*
            sample->pose.v[0] += (drand48() * 4 * data->dispersion.v[0]) - (2 * data->dispersion.v[0]);
            sample->pose.v[1] += (drand48() * 4 * data->dispersion.v[1]) - (2 * data->dispersion.v[1]);
            sample->pose.v[2] += (drand48() * 6 * data->dispersion.v[2]) - (3 * data->dispersion.v[2]);
            //sample->pose.v[2] =  pf_ran_gaussian(3.14);
            //*/
            /*
            sample->pose.v[0] += (drand48() *4 ) - 2;
            sample->pose.v[1] += (drand48() *4 ) - 2;
            sample->pose.v[2] += (drand48() *2) -1 ;
            */

            sample->weight = 1.0 / set->sample_count;
        }
    }


}

void AmclgnssSensor::er_gr(double kld, pf_t *pf, gnssSensorData *gnss_data, double beta,  analysis_t *data){
    int reset_count = 3;
    pf_sample_set_t *set;
    pf_sample_t *sample;
    int particle_num;
    pf_vector_t pose;
    const double gnss_sigma = 2.0;
    set = pf->sets + pf->current_set;
    int gr_sample_count = beta * set->sample_count;
    if(kld > 20 && (data->dispersion.v[0]+data->dispersion.v[0]) > 10){
        //GRリセットを行う
        //std::cout << "GR !!!!!!!!!!!!" << std::endl;
        for(int i=0; i<gr_sample_count; i++){
            particle_num = rand() % set->sample_count + 1;
            sample = set->samples + particle_num;
            pose = sample->pose;

            sample->pose.v[0] = gnss_data->x + pf_ran_gaussian(gnss_sigma);
            sample->pose.v[1] = gnss_data->y + pf_ran_gaussian(gnss_sigma);
            sample->pose.v[2] = pf_ran_gaussian(3.14);
        }
    }else{
        //膨張リセットを行う
        //std::cout << "ER !!!!!!!!!!!!" << std::endl;
        double reset_limit = ((int)data->dispersion.v[0] + (int)data->dispersion.v[1]) / 2;
        if(reset_count >= reset_limit){
            for(int i=0; i<set->sample_count; i++ ){
                sample = set->samples + i;
                sample->pose.v[0] += (drand48() * 4 * data->dispersion.v[0]) - (2 * data->dispersion.v[0]);
                sample->pose.v[1] += (drand48() * 4 * data->dispersion.v[1]) - (2 * data->dispersion.v[1]);
                sample->pose.v[2] += (drand48() * 2 * data->dispersion.v[2]) - (1 * data->dispersion.v[2]);
                sample->weight = 1.0 / set->sample_count;
            }
        }
    }

}
