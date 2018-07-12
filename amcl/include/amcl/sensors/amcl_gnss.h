#ifndef gnss_SENSOR_H
#define gnss_SENSOR_H


#include "amcl_sensor.h"

namespace amcl{

class gnssSensorData : public AMCLSensorData{
public:
    double x;
    double y;
    int position_covariance_type;
};

class AmclgnssSensor : public gnssSensorData{
private:
    static double normal_distribution(double diff, double dispersion);
    static double gnssModel(AmclgnssSensor *gnss_data, pf_sample_set_t* set);
public:
    AmclgnssSensor();
    //~AmclgnssSensor(){};

    virtual bool gnssSensorUpdata(pf_t *pf, gnssSensorData *gnss_data);
    virtual bool gnssSensor_reseting(pf_t *pf, gnssSensorData *gnss_data, double sample_num);
    virtual double gnssPfKLD(pf_vector_t pf, gnssSensorData *gnss_data);
    virtual void er(pf_t *pf, analysis_t *data);
    virtual void er_gr(double kld, pf_t *pf, gnssSensorData *gnss_data, double beta,  analysis_t *data);
};

}
#endif
