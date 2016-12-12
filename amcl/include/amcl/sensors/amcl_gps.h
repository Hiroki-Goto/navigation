#ifndef GPS_SENSOR_H
#define GPS_SENSOR_H


#include "amcl_sensor.h"

namespace amcl{

class GPSSensorData : public AMCLSensorData{
public:
    double x;
    double y;
    int position_covariance_type;
};

class AmclGPSSensor : public GPSSensorData{
private:
    static double normal_distribution(double diff, double dispersion);
    static double GPSModel(AmclGPSSensor *gps_data, pf_sample_set_t* set);
public:
    AmclGPSSensor();
    //~AmclGPSSensor(){};

    virtual bool gpsSensorUpdata(pf_t *pf, GPSSensorData *gps_data);
};

}
#endif
