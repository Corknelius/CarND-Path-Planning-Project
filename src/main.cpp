#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("}");
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for(int i = 0; i < maps_x.size(); i++)
    {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x,y,map_x,map_y);
        if(dist < closestLen)
        {
            closestLen = dist;
            closestWaypoint = i;
        }

    }

    return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

    int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double heading = atan2((map_y-y),(map_x-x));

    double angle = fabs(theta-heading);
    angle = min(2*pi() - angle, angle);

    if(angle > pi()/4)
    {
        closestWaypoint++;
        if (closestWaypoint == maps_x.size())
        {
            closestWaypoint = 0;
        }
    }

    return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

    int prev_wp;
    prev_wp = next_wp-1;
    if(next_wp == 0)
    {
        prev_wp  = maps_x.size()-1;
    }

    double n_x = maps_x[next_wp]-maps_x[prev_wp];
    double n_y = maps_y[next_wp]-maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
    double proj_x = proj_norm*n_x;
    double proj_y = proj_norm*n_y;

    double frenet_d = distance(x_x,x_y,proj_x,proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000-maps_x[prev_wp];
    double center_y = 2000-maps_y[prev_wp];
    double centerToPos = distance(center_x,center_y,x_x,x_y);
    double centerToRef = distance(center_x,center_y,proj_x,proj_y);

    if(centerToPos <= centerToRef)
    {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for(int i = 0; i < prev_wp; i++)
    {
        frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
    }

    frenet_s += distance(0,0,proj_x,proj_y);

    return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int prev_wp = -1;

    while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
    {
        prev_wp++;
    }

    int wp2 = (prev_wp+1)%maps_x.size();

    double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s-maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
    double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

    double perp_heading = heading-pi()/2;

    double x = seg_x + d*cos(perp_heading);
    double y = seg_y + d*sin(perp_heading);

    return {x,y};

}

int find_lane(double car_d) {
    // Left Lane: 0
    // Center Lane: 1
    // Right Lane: 2
    // Offroad: -1
    int lane = 0;
    if (car_d >= 0.0 && car_d <= 12.0) {
        lane = car_d/4;
    }
    else {
        lane = -1;
    }
    return lane;
}

int main() {
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;

    ifstream in_map_(map_file_.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line)) {
        istringstream iss(line);
        double x;
        double y;
        float s;
        // normal components to waypoints
        float d_x;
        float d_y;
        iss >> x;
        iss >> y;
        iss >> s;
        iss >> d_x;
        iss >> d_y;
        map_waypoints_x.push_back(x);
        map_waypoints_y.push_back(y);
        map_waypoints_s.push_back(s);
        map_waypoints_dx.push_back(d_x);
        map_waypoints_dy.push_back(d_y);
    }

    // TODO: have some vars here
    double target_v = 0.;
    double target_d = 6.0;

    h.onMessage([&target_v, &target_d, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](
            uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
            uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        //auto sdata = string(data).substr(0, length);
        //cout << sdata << endl;
        if (length && length > 2 && data[0] == '4' && data[1] == '2') {

            auto s = hasData(data);

            if (s != "") {
                auto j = json::parse(s);

                string event = j[0].get<string>();

                if (event == "telemetry") {
                    // j[1] is the data JSON object

                    // Main car's localization Data
                    double car_x = j[1]["x"];
                    double car_y = j[1]["y"];
                    double car_s = j[1]["s"];
                    double car_d = j[1]["d"];
                    double car_yaw = j[1]["yaw"];
                    double car_speed = j[1]["speed"];
                    int car_lane = find_lane(car_d);

                    // Previous path data given to the Planner
                    auto previous_path_x = j[1]["previous_path_x"];
                    auto previous_path_y = j[1]["previous_path_y"];
                    // Previous path's end s and d values
                    double end_path_s = j[1]["end_path_s"];
                    double end_path_d = j[1]["end_path_d"];

                    // Sensor Fusion Data, a list of all other cars on the same side of the road.
                    auto sensor_fusion = j[1]["sensor_fusion"];

                    json msgJson;

                    vector<double> next_x_vals;
                    vector<double> next_y_vals;


                    //------------------
                    // Begin Edits here:
                    //------------------

                    int path_size = 40;

                    if (previous_path_x.size() > 0)
                    {
                        car_s = end_path_s;
                    }

                    // Step 1: Integrate Sensor Fusion
                    //Initialize vectors to store traffic info
                    vector<double> traffic_same_s;
                    vector<double> traffic_same_v;
                    
                    vector<double> traffic_left_s;
                    vector<double> traffic_left_v;
                    
                    vector<double> traffic_right_s;
                    vector<double> traffic_right_v;
                    
                    double closest_same_front_dist = 9999.99;
                    double closest_same_back_dist = 9999.99;
                    double closest_same_front_v = 0.0;
                    
                    double closest_left_front_dist = 9999.99;
                    double closest_left_back_dist = 9999.99;
                    double closest_left_front_v = 0.0;
                    double closest_left_back_v = 0.0;
                    
                    double closest_right_front_dist = 9999.99;
                    double closest_right_back_dist = 9999.99;
                    double closest_right_front_v = 0.0;
                    double closest_right_back_v = 0.0;
                    
                    // Detect other cars
                    for (int i = 0; i < sensor_fusion.size(); i++) {
                        double traffic_id = sensor_fusion[i][0];
                        double traffic_x = sensor_fusion[i][1];
                        double traffic_y = sensor_fusion[i][2];
                        double traffic_vx = sensor_fusion[i][3];
                        double traffic_vy = sensor_fusion[i][4];
                        double traffic_s = sensor_fusion[i][5];
                        double traffic_d = sensor_fusion[i][6];
                        double traffic_v = sqrt(traffic_vx * traffic_vx + traffic_vy * traffic_vy);
                        double traffic_distance = sensor_fusion[i][6];
                        int traffic_lane = find_lane(traffic_d);
                        
                        if (traffic_distance < 300) {
                            // Sort traffic into their respective lanes
                            if (car_lane == traffic_lane) {
                                // Same Lanes
                                traffic_same_s.push_back(traffic_s);
                                traffic_same_v.push_back(traffic_v);
                            }
                            //else if ((traffic_lane >= (car_lane + 1)) and traffic_lane < 2) {
                            else if ((traffic_lane == car_lane + 1) and traffic_lane <= 2) {
                                // Cars to the Right
                                traffic_right_s.push_back(traffic_s);
                                traffic_right_v.push_back(traffic_v);
                            }
                            //else if ((traffic_lane <= (car_lane - 1)) and traffic_lane > 0) {
                            else if ((traffic_lane == (car_lane - 1)) and traffic_lane >= 0) {
                                // Cars to the Left
                                traffic_left_s.push_back(traffic_s);
                                traffic_left_v.push_back(traffic_v);
                            }
                        }
                    }
                    
                    cout<<endl;
                    cout<<"Traffic_mylane:"<< traffic_same_s.size()<<endl;
                    cout<<"Traffic_myleft:"<< traffic_left_s.size()<<endl;
                    cout<<"Traffic_myright:"<< traffic_right_s.size()<<endl;;
                    cout<<endl;
                    
                    
                    
                    // Find closest vehicles in each lane
                    if (traffic_same_s.size() > 0) {
                        for (int i=0; i < traffic_same_s.size(); i++) {
                            double distance_same = car_s - traffic_same_s[i];
                            if (distance_same > 0 and abs(distance_same) < abs(closest_same_back_dist)) {
                                // Cars behind
                                closest_same_back_dist = abs(distance_same);
                            }
                            else if (distance_same < 0 and abs(distance_same) < abs(closest_same_back_dist)) {
                                // Cars in front
                                closest_same_front_dist = abs(distance_same);
                                closest_same_front_v = traffic_same_v[i];
                            }
                        }
                    }
                    
                    if (traffic_left_s.size() > 0) {
                        for (int i=0; i < traffic_left_s.size(); i++) {
                            double distance_left = car_s - traffic_left_s[i];
                            if (distance_left > 0 and abs(distance_left) < abs(closest_left_back_dist)) {
                                // Cars behind
                                closest_left_back_dist = abs(distance_left);
                                closest_left_back_v = traffic_left_v[i];
                            }
                            else if (distance_left < 0 and abs(distance_left) < abs(closest_left_back_dist)) {
                                // Cars in front
                                closest_left_front_dist = abs(distance_left);
                                closest_left_front_v = traffic_left_v[i];
                            }
                        }
                    }
                    
                    if (traffic_right_s.size() > 0) {
                        for (int i=0; i < traffic_right_s.size(); i++) {
                            double distance_right = car_s - traffic_right_s[i];
                            if (distance_right > 0 and abs(distance_right) < abs(closest_right_back_dist)) {
                                // Cars behind
                                closest_right_back_dist = abs(distance_right);
                                closest_right_back_v = traffic_right_v[i];
                            }
                            else if (distance_right < 0 and abs(distance_right) < abs(closest_right_back_dist)) {
                                // Cars in front
                                closest_right_front_dist = abs(distance_right);
                                closest_right_front_v = traffic_right_v[i];
                            }
                        }
                    }
                    
                    /*
                    cout<<endl;
                    cout<<"Closest Same Front:"<< closest_same_front_dist<<endl;
                    cout<<"Closest Same Front Vel:"<< closest_same_front_v<<endl;
                    cout<<"Closest Same Back:"<< closest_same_back_dist<<endl;

                    cout<<"Closest Left Front:"<< closest_left_front_dist<<endl;
                    cout<<"Closest Left Front Vel:"<< closest_left_front_v<<endl;
                    cout<<"Closest Left Back:"<< closest_left_back_dist<<endl;
                    cout<<"Closest Left Back Vel:"<< closest_left_back_v<<endl;
                    
                    cout<<"Closest Right Front:"<< closest_right_front_dist<<endl;
                    cout<<"Closest Right Front Vel:"<< closest_right_front_v<<endl;
                    cout<<"Closest Right Back:"<< closest_right_back_dist<<endl;
                    cout<<"Closest Right Back Vel:"<< closest_right_back_v<<endl;
                    cout<<endl;
                    */
                    
                    // Step 2: Calculate Costs
                    // Initialize costs: 0 = keep, 1 = left, 2 = right
                    vector<int> costs = {0,0,0};
                    
                    double SPEED_LIMIT = 43.8;  // Not exactly 43.8, need to figure out why it's closer to 47.
                    double FREE_ACCEL = 0.2;    // When no cars in front
                    double CLOSE_ACCEL =0.05;   // When other cars are nearby
                    
                    // Make going out of bounds impossible
                    if (car_lane == 0) {
                        costs[1] = 9999;
                    } else if (car_lane == 2) {
                        costs[2] = 9999;
                    }
                    // Calc cost for staying in lane
                    costs[0] -= (closest_same_front_dist/5) - 8;

                    // Calc cost for changing left
                    if (closest_left_front_dist > 40 and closest_left_back_dist > 40) {
                        costs[1] -= 5; // Incentivise going into an open lane, but not as much as keeping lane
                    } else if ((closest_left_back_dist < 20) or (closest_left_front_dist < 20)) {
                        costs[1] += 10; //Penalize changing lanes with no room
                    }
                    if (closest_left_front_v > closest_same_front_v) {
                        costs[1] -= 3; // Incentivise going into a faster lane
                    }
                    
                    // Calc cost for changing right
                    if (closest_right_front_dist > 40 and closest_right_back_dist > 40) {
                        costs[2] -= 5; // Incentivise going into an open lane, but not as much as keeping lane
                    } else if ((closest_right_back_dist < 20) or (closest_right_front_dist < 20)) {
                        costs[2] += 10; //Penalize changing lanes with no room
                    }
                    if (closest_right_front_v > closest_same_front_v) {
                        costs[2] -= 3; // Incentivise going into a faster lane
                    }
                    
                    // Step 3: Process Costs into Behaviors
                    // Pick minimum cost to be the behavior
                    int behavior = 0;
                    int min_cost = 9999;
                    
                    for (int i = 0; i < 3; i++) {
                        int temp_cost = costs[i];
                        if (temp_cost < min_cost) {
                            behavior = i;
                            min_cost = costs[i];
                        }
                    }
                    
                    cout<<endl;
                    cout<<"Cost Keep:"<< costs[0] <<endl;
                    cout<<"Cost Change Left:"<< costs[1] <<endl;
                    cout<<"Cost Change Right:"<< costs[2] <<endl;
                    cout<<"Behavior: "<<  behavior <<endl;
                    cout<<"Min Cost: "<<  min_cost <<endl;
                    cout<<endl;
                    
                   
                    if (behavior == 0) {
                        // Set Speed
                        if (closest_same_front_dist < 30) {
                            cout << "Closest front dist: " << closest_same_front_dist << endl;
                            cout << "Closest front Vel: " << closest_same_front_v << endl;
                            if (car_speed >= closest_same_front_v) {
                                target_v -= CLOSE_ACCEL;
                            } else if (car_speed < closest_same_front_v) {
                                if (car_speed < SPEED_LIMIT) {
                                    target_v += CLOSE_ACCEL;
                                } else {
                                    cout << "Matching closest front car's speed." << endl;
                                }
                            }
                        } else if (30 < closest_same_front_dist < 75) {
                            if (car_speed < SPEED_LIMIT) {
                                target_v += CLOSE_ACCEL;
                            }
                        } else {
                            if (car_speed < SPEED_LIMIT) {
                                target_v += FREE_ACCEL;
                            } else {
                                cout << "Maintaining speed: " << car_speed << endl;
                            }
                        }

                        target_d = 2.0 + (car_lane * 4.0);
                        
                    } else if (behavior == 1) {
                        // Set Position
                        if (car_lane == 0) {
                            target_d = 2.0; // Shouldn't be possible, but just in case...
                        } else if (car_lane == 1) {
                            target_d = 2.0;
                        } else if (car_lane == 2) {
                            target_d = 6.0;
                        } else {
                            target_d = car_d;
                        }
                    } else if (behavior == 2) {
                        // Set Position
                        if (car_lane == 0) {
                            target_d = 6.0;
                        } else if (car_lane == 1) {
                            target_d = 10.0;
                        } else if (car_lane == 2) {
                            target_d = 6.0;  // Shouldn't be possible, but just in case...
                        } else {
                            target_d = car_d;
                        }
                    } else {
                        cout<<"Something is wrong, behavior is : " << behavior;
                    }

                    // Step 4: Generate path points
                    
                    vector<double> spline_x;
                    vector<double> spline_y;

                    double ref_x = car_x;
                    double ref_y = car_y;
                    double ref_yaw = deg2rad(car_yaw);

                    if(previous_path_x.size() < 2)
                    {

                        double prev_car_x = car_x - cos(car_yaw);
                        double prev_car_y = car_y - sin(car_yaw);

                        spline_x.push_back(prev_car_x);
                        spline_x.push_back(car_x);

                        spline_y.push_back(prev_car_y);
                        spline_y.push_back(car_y);

                    }
                    else
                    {
                        ref_x = previous_path_x[previous_path_x.size() - 1];
                        ref_y = previous_path_y[previous_path_x.size() - 1];

                        double ref_x_prev = previous_path_x[previous_path_x.size() - 2];
                        double ref_y_prev = previous_path_y[previous_path_x.size() - 2];
                        ref_yaw = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);
                        spline_x.push_back(ref_x_prev);
                        spline_x.push_back(ref_x);

                        spline_y.push_back(ref_y_prev);
                        spline_y.push_back(ref_y);

                    };

                    /*
                     * Generate three spaced out points and draw a spline through them
                     */

                    vector<double> next_wp0 = getXY(car_s + 30, target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    vector<double> next_wp1 = getXY(car_s + 60, target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    vector<double> next_wp2 = getXY(car_s + 90, target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

                    spline_x.push_back(next_wp0[0]);
                    spline_x.push_back(next_wp1[0]);
                    spline_x.push_back(next_wp2[0]);

                    spline_y.push_back(next_wp0[1]);
                    spline_y.push_back(next_wp1[1]);
                    spline_y.push_back(next_wp2[1]);


                    /*
                     * Convert points to frame of ref relative to car
                     */

                    for (int i = 0; i < spline_x.size(); i++) {
                        double shift_x = spline_x[i] - ref_x;
                        double shift_y = spline_y[i] - ref_y;

                        spline_x[i] = (shift_x*cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
                        spline_y[i] = (shift_x*sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
                    }

                    tk::spline s;

                    s.set_points(spline_x, spline_y);

                    // Add unprocessed points from previous path to current path

                    for(int i = 0; i < previous_path_x.size(); i++)
                    {
                        next_x_vals.push_back(previous_path_x[i]);
                        next_y_vals.push_back(previous_path_y[i]);
                    }


                    // Set horizon
                    double target_x = 30.;
                    double target_y = s(target_x);
                    double target_dist = sqrt((target_x)*(target_x) + (target_y)*(target_y));

                    double x_add_on = 0;

                    for (int i = 0; i < path_size - previous_path_x.size(); i++) {

                        double N = (target_dist/(.02*target_v));
                        double x_point = x_add_on + (target_x)/N;
                        double y_point = s(x_point);

                        x_add_on = x_point;

                        // Convert coords from car's frame of ref

                        double x_ref = x_point;
                        double y_ref = y_point;

                        x_point = (x_ref * cos(ref_yaw)-y_ref*sin(ref_yaw));
                        y_point = (x_ref * sin(ref_yaw)+y_ref*cos(ref_yaw));

                        x_point += ref_x;
                        y_point += ref_y;

                        next_x_vals.push_back(x_point);
                        next_y_vals.push_back(y_point);

                    }


                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;

                    auto msg = "42[\"control\"," + msgJson.dump() + "]";

                    //this_thread::sleep_for(chrono::milliseconds(1000));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

                }
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
