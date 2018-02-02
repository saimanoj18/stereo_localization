#include <iostream>
#include <ros/ros.h>
#include "CamLocalization.h"

using namespace std;



int main(int argc, char **argv)
{
    ros::init(argc, argv, "camera_localization");
    ros::start();
    CamLocalization camloc;
    
    ros::Rate r(30);
    while (ros::ok())
    {
        camloc.Refresh();
        ros::spinOnce();
        r.sleep();
    }
    ros::shutdown();
    return 0;
}
