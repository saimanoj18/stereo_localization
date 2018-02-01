#include <iostream>
#include <ros/ros.h>
#include "CamLocalization.h"

using namespace std;

int64_t
timestamp_now (void)
{
    struct timeval tv;
    gettimeofday (&tv, NULL);
    return (int64_t) tv.tv_sec * 1000000 + tv.tv_usec;
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "camera_localization");
    ros::start();
    CamLocalization camloc;
    
    //for time computation
    int64_t start_time;
    int64_t end_time;
    
    ros::Rate r(30);
    while (ros::ok())
    {
        start_time = timestamp_now ();
        camloc.Refresh();
        end_time = timestamp_now ();
        float time_diff = (end_time - start_time)/1000000.0;
        if(time_diff>0.001){
        camloc.write_times("Elapsed_times.txt", time_diff);
        cout<<"Elapsed time: %"<<time_diff<<" secs\n"<<endl;
        }
        ros::spinOnce();
        r.sleep();
    }
    ros::shutdown();
    return 0;
}
