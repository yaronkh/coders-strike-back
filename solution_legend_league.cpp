#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <memory>
#include <utility>
#include <tuple>
#include <bits/stdc++.h>

using namespace std;

#define PI 3.14159265

double radians(double ang) {
   return ang * PI / 180.;
}

double degrees(double rad) {
   return rad * 180.0 / PI;
}

template <typename T=int>
T clip(T val, T v1, T v2) {
      if (val < v1)
         val = v1;
      else if(val > v2)
         val = v2;
      return val;
}

template <typename T>
struct Coord {
   T x;
   T y;

   friend ostream & operator << (ostream &out, const Coord &c) {
      out << "(" << c.x << "," << c.y << ")";
      return out;
   }

   friend Coord operator - (const Coord c0,const Coord &c1) {
      return {c0.x - c1.x, c0.y - c1.y};
   }

   double norm(void) const {
      return sqrt(x*x + y*y);
   }

   double dist_pnts(const Coord &c) {
      return (*this - c).norm();
   }

   T norm2(void) const {
      return x*x + y*y;
   }

   int ang(void) {
      return atan2(y, x) * 180./PI;
   }

   static T det(const Coord &c1, const Coord &c2) {
      return c1.x * c2.y - c1.y * c2.x;
   }

   static T dot(const Coord &c1, const Coord &c2) {
      return c1.x * c2.x + c1.y * c2.y;
   }

   T dot(const Coord &c) const {
      return dot(*this, c);
   }

   Coord operator / (T t) const {
      return {x / t, y / t};
   }

   Coord operator * (T t) {
      return {x * t, y * t};
   }

   friend Coord operator + (const Coord &p1, const Coord &p2) {
      return {p1.x + p2.x, p1.y + p2.y};
   }

   T cross(const Coord &c) {
      return x * c.y - y * c.x;
   }

   Coord unit_vec(void) const {
      return (*this)/norm();
   }

   double angle_between(const Coord &v2) {
      auto v1_u = unit_vec();
      auto v2_u = v2.unit_vec();
      auto c = clip(dot(*this, v2), -1., 1.);
      auto ret = degrees(acos(c));
      auto cr = v1_u.cross(v2_u);
      return (cr >= 0) ? ret : -ret;
   }
};

typedef Coord<int>   icoord;
typedef Coord<double> dcoord;

dcoord to_dcoord(const icoord &d) {
   return {double(d.x), double(d.y)};
}

icoord to_icoord(const dcoord &d) {
   return {int(round(d.x)), int(round(d.y))};
}

template <typename T=int>
double relAngle(const Coord<T> &p1, const Coord<T> &p2) {
   return atan2(Coord<T>::det(p1, p2), Coord<T>::dot(p1, p2));
}

dcoord unit_coord(double angle)
{
   auto ang_rad = radians(angle);
   return {cos(ang_rad), sin(ang_rad)};
}

template <typename T = int>
struct Line {
   Coord<T> p1;
   Coord<T> p2;

   T point_to_line_dist(Coord<T> p) {
      T x1 = p1.x;
      T y1 = p1.y;
      T x2 = p2.x;
      T y2 = p2.y;
      T x0 = p.x;
      T y0 = p.y;
      return -((y2 - y1)*x0 - (x2 - x1) *y0 +x2*y1 - y2*x1)/sqrt((y2 - y1)*(y2 - y1) + (x2 - x1)*(x2 - x1));
   }

};

typedef Line<double>   fline;

void line_ABC(const dcoord &p1, const dcoord &p2, double &A, double &B, double &C)
{
   A = p2.y - p1.y;
   B = p1.x - p2.x;
   C = A * p1.x + B * p1.y;
}

dcoord line_intersect(const dcoord &p1, const dcoord &p2, const dcoord &q1, const dcoord &q2)
{
   double A1, B1, C1, A2, B2, C2;
   line_ABC(p1, p2, A1, B1, C1);
   line_ABC(q1, q2, A2, B2, C2);
   double det = A1 * B2 - A2 * B1;
   if (det == 0.0)
      throw "attempt to find intersect of parallel lines";
   return {(B2 * C1 - B1 * C2) / det, (A1 * C2 - A2 * C1) / det};
}

void test_line_intersect()
{
   dcoord p1 = {0., 5.0}, p2 = {5.0, 0.0};
   dcoord q1 = {0.0, 0.0}, q2 = {5.0, 5.0};
   cerr << "line intersect:" << p1 << p2 << " and " << q1 << q2 << endl;
   cerr << "result:" << line_intersect(p1, p2, q1, q2) << endl;
}

dcoord rotate(const dcoord &v, double tet_degrees, const dcoord &o = {0, 0})
{
   auto vx = v.x - o.x;
   auto vy = v.y - o.y;
   auto teta = radians(tet_degrees);
   auto costet = cos(teta);
   auto sintet = sin(teta);
   auto v10 = vx * costet;
   auto v11 = vx * sintet;
   auto v20 = vy * sintet;
   auto v21 = vy * costet;
   return {v10 - v20 + o.x, v11 + v21 + o.y};
}

void test_rotate()
{
   double teta = 36.8698;
   dcoord B = {6., 0.};
   dcoord A = {2., 3.};
   dcoord C = rotate(B, teta, A);
   cerr << "rotation of " << B << " around " << A << " gives " << C << endl;
}

struct fcirc {
   dcoord c;
   double rad;

   friend ostream & operator << ( ostream &out, const fcirc &cr) {
      out << "center: " << cr.c << "rad: " << cr.rad;
      return out;
   }
};

fcirc find_rad_from_two_points_and_tangent(const dcoord &p1, const dcoord &tang, const dcoord &p2)
{
    /*
    p1 - point on the circle
    tang - direction vector relative to p1
    p2 - point on the circle
   */
   fcirc res;
   auto pepend1 = rotate(tang, 90.0) + p1;
   auto p3 = (p1 + p2) / 2.0;
   auto pepend2 = rotate(p2, 90, p3);
   res.c = line_intersect(p1, pepend1, p3, pepend2);
   res.rad = (res.c - p1).norm();
   return res;
}

double angleabs2angle(double angle_abs, const dcoord &target, const dcoord &pos)
{
   auto angle_uv = unit_coord(radians(angle_abs));
   auto t = target - pos;
   return angle_uv.angle_between(t);
}

void test_angleabs2angle()
{
   dcoord pos = {10.0, 10.0};
   dcoord target = {500., 500.};
   double angle_abs = 45.0;
   cerr << "pos=" << pos << " target=" << target << " face=" << angle_abs << "ang_rel=" <<
      angleabs2angle(angle_abs, target, pos) << endl;
}

/*
 * In [3]: s.find_rad_from_two_points_and_tangent((0.0, 5.0), (1., 1.), (5.0, 0.0))
Out[3]: (3.5355339059327378, (2.4999999999999996, 2.4999999999999996))

 */
void test_circ_from_2points_and_rad()
{
   dcoord p1 = {0.0, 5.0};
   dcoord tang = {1., 1.};
   dcoord p2 = {5.0, 5.0};
   cerr << "circ from " << p1 << " tangent " << tang << " and " << p2 << endl;
   cerr << "result: " << find_rad_from_two_points_and_tangent(p1, tang, p2) << endl;
}

void test_angle_between()
{
   dcoord v1 = {10.0, 0};
   dcoord v2 = {0., 11};
   cerr << "angle v1 to v2 "<< v1.angle_between(v2) << endl;
   cerr << "angle v2 to v1 "<< v2.angle_between(v1) << endl;
}

double location_along_segment(const dcoord &p, const dcoord &q, const dcoord &x)
{
   auto v = q - p;
   auto s = x - p;
   return s.dot(v) / v.dot(v);
}

/*In [2]: s.location_along_segment((0.0, 0.0), (5.0, 5.0), (-1.0, 8.0))
Out[2]: 0.7
*/
void test_location_along_segment()
{
   dcoord p = {0., 0.};
   dcoord q = {5.0, 5.0};
   dcoord x = {-1.0, 8.0};
   cerr << "seg " << p << q << " point " << x << " loc along seg=" << location_along_segment(p,q,x) << endl;
}

typedef icoord ivec;

class Planner;

struct PodParams {
   static constexpr int max_thrust = 200;
   static constexpr int rot_vel = 18;
   static constexpr double fric_thruts_to_thrust_ratio = 0.85;

   unique_ptr<Planner> planner = make_unique<Planner>();

};
struct ArenaParams {
   static constexpr int station_rad = 600;
   static constexpr double friction_fac = 0.85;
   static constexpr int station_tolerance = 450;

   bool start_with_boost = false;
   int defender_dist_spare = 1000;
   double gtKp = 1.0;
   double gtKi = 0.;
   double gtKd = 0;
};

ArenaParams HostileParams() {
   return ArenaParams();
};

ArenaParams PyramidParams() {
   return ArenaParams();
};

ArenaParams TriangleParams() {
   return ArenaParams();
};

ArenaParams DaltonParams() {
   auto ret = ArenaParams();
   ret.gtKp = 0.5;
   ret.start_with_boost = true;
   return ret;
};

ArenaParams ArrowParams() {
   return ArenaParams();
};

ArenaParams MacbilitParams() {
   return ArenaParams();
};

ArenaParams ShoshParams() {
   return ArenaParams();
};

ArenaParams TilParams() {
   return ArenaParams();
};

ArenaParams TrapezParams() {
   return ArenaParams();
};

ArenaParams MehumashParams() {
   return ArenaParams();
};

ArenaParams TrampolineParams() {
   return ArenaParams();
};

ArenaParams ZigzagParams() {
   return ArenaParams();
};

typedef tuple<icoord, int, bool, bool> instruction;

class Runner;

class Planner {
public:

   void reset(Runner *runner_) { runner = runner_; }

   virtual void plan(icoord target, Runner *runner_) {
      runner = runner_;
   }

   virtual instruction act(void) {
      icoord target = {0, 0};
      return make_tuple(target, 0, false, false);
   }

protected:
   Runner *runner = nullptr;
};

class Predictor
{
public:
    Predictor(Runner *runner);

private:


};

class Runner;

class GoToTargetRegulator: public Planner {
public:
   void reset(double tolerance_, Runner *runner_) {
      tolerance = tolerance_;
      runner = runner_;
      error = 0.0;
      ie = 0.0;
      dedt = 0.0;
      last_e = 0.0;
      ttarget = {0, 0};
      thrust_reducing = PodParams::max_thrust;
      is_pointing = false;
      is_reducing = false;
   }
   virtual instruction act(const dcoord &itarget);
   void calc_target_vel(const dcoord &target, double &vt, double &alpha);
   bool is_pointing_target(const dcoord &target);
   int regulate_thrust(const dcoord &target, double angle);
   instruction reduce(const dcoord &target, int thrust);

   double error = 0.;
   double ie = 0.0;
   double last_e = 0.0;
   double dedt  = 0.0;
   Runner *runner;
   double tolerance = 0.0;
   bool is_pointing = false;
   int thrust_reducing = PodParams::max_thrust;
   bool is_reducing = false;
   dcoord ttarget = {0, 0};
   int thrust = PodParams::max_thrust;
};
/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/

class Runner
{
public:
   Runner(int num_laps, bool is_me) {


   }

   dcoord direction(void) {
      auto a_vel = abs_vel();
      if (a_vel == 0.0)
         a_vel = 1.0;
      return {cvel().x/a_vel, cvel().y / a_vel};
   }

   static void push_me(unique_ptr<Runner> r) {
       me.push_back(std::move(r));
   }

   static void push_opponent(unique_ptr<Runner> r) {
       other.push_back(std::move(r));
   }

   icoord station(int n) {
      auto itr = stations.begin();
      for (int i = 0; i < n; ++i)
         ++itr;
      return *itr;
   }

   icoord &cpos(void) {
      return pos.back();
   }

   int cangle(void) {
      return angles.back();
   }

   double abs_vel(void) {
      return vels.back().norm();
   }

   ivec &cvel(void) {
      return vels.back();
   }

   double pod_deflection(void) {
      double al = radians(cangle());
      dcoord face = {cos(al), sin(al)};
      return to_dcoord(cvel()).angle_between(face);
   }

public:
   static                  vector<unique_ptr<Runner> > me;
   static                  vector<unique_ptr<Runner> > other;

   int                     num_laps = 0;
   ArenaParams             arena_params;
   PodParams               pod_params;
   list<icoord>            pos;
   list<int>               angles;
   list<ivec>              vels;
   list<int>               thrusts;
   int                     next_cp = 0;
   int                     prev_cp = 0;
   list<icoord>            stations;
   int                     time_left_to_punch = 0;
   int                     passed_stations = -1;
   bool                    is_me = false;
   bool                    boost = false;
   int                     boost_turns = 0;
   shared_ptr<Predictor>   predictor;
   list<icoord>            predictions;
};

instruction GoToTargetRegulator::reduce(const dcoord &target, int thrust)
{
      if (!is_reducing) {
         is_reducing =  true;
         thrust_reducing = thrust - 30;
      } else {
         thrust_reducing -= 30;
      }
      if (thrust_reducing < 80)
         thrust_reducing = 80;
      return make_tuple(to_icoord(target), thrust_reducing, false, false);
}

instruction GoToTargetRegulator::act( const dcoord &target)
{
   if (is_pointing_target(target)) {
      is_pointing = true;
      auto ret = make_tuple(to_icoord(ttarget), PodParams::max_thrust, false, false);
      cerr << "POINTED TO TARGET" << endl;
      return ret;
   }
   is_pointing = false;
   auto dir_ = runner->direction();
   auto pt = target - to_dcoord(runner->cpos());
   thrust = PodParams::max_thrust;
   auto v = runner->abs_vel();
   auto angle = runner->cangle();
   if (v > 1.0)
      thrust = regulate_thrust(target, angle);
   error = 0.0;
   if (dir_.x != 0.0 || dir_.y != 0.0)
      error = pt.angle_between(dir_);
   if (fabs(error) > 89.) {
      reset(tolerance, runner);
      return reduce(target, thrust);
   }
   ie += error;
   dedt = last_e - error;
   auto kp = runner->arena_params.gtKp;
   auto ki = runner->arena_params.gtKi;
   auto kd = runner->arena_params.gtKd;
   auto c_angle = -(kp*error + ki*ie + kd*dedt);
   last_e = error;
   auto pt1 = rotate(pt, c_angle);
   auto a_pt = rotate(pt, angle);
   auto diff = fabs(pt1.angle_between(a_pt));
   dcoord res = to_dcoord(runner->cpos()) + pt1;
   double d = to_dcoord(runner->cpos()).dist_pnts(target);
   cerr << "DDD=" << d << endl;
   if (diff > 18.0 && d < 3500) {
      return reduce(res, thrust);
   }
   is_reducing = false;
   return make_tuple(to_icoord(res), thrust, false, false);
}

void GoToTargetRegulator::calc_target_vel(const dcoord &target, double &vt, double &alpha)
{
   auto p1 = to_dcoord(runner->cpos());
   auto tg = to_dcoord(runner->cvel()).unit_vec();
   auto circ = find_rad_from_two_points_and_tangent(p1, tg, target);
   auto rad = circ.rad + 600.;
   alpha = fabs(radians(runner->pod_deflection()));
   vt = rad * (1.0 - ArenaParams::friction_fac) * tan(alpha);
   static const double omega = radians(PodParams::rot_vel);
   double vt2 = omega * rad;
   if (vt2 < vt) {
      cerr << "reducing velocity to meet max radial vel" << endl;
      vt = vt2;
   }
}

bool GoToTargetRegulator::is_pointing_target(const dcoord &target)
{
   if (runner->abs_vel() < 2.0) {
      cerr << "VELOCITY TOO SMALL TO TELL POINTING" << endl;
      return false;
   }
   auto direc = to_dcoord(runner->cvel()).unit_vec();
   auto p1 = to_dcoord(runner->cpos());
   auto p2 = p1 + direc;
   auto perp = rotate(direc, 90);
   const dcoord &p3 = target;
   auto p4 = p3 + perp;
   ttarget = line_intersect(p1, p2, p3, p4);
   double fac = location_along_segment(p1, p2, ttarget);
   double d = ttarget.dist_pnts(target);
   return d <= tolerance && fac >= 0.;
}

int GoToTargetRegulator::regulate_thrust(const dcoord &target, double angle)
{
   double pod_def = runner->pod_deflection();
   double v_tar = 0;
   double alpha = 0.;
   if (angle < 1.0 and fabs(pod_def) < 1.0) {
      v_tar = INT_MAX;
      alpha = 0.0;
   } else {
      calc_target_vel(target, v_tar, alpha);
   }
   double thrust_ = v_tar * (1 - ArenaParams::friction_fac) / cos(alpha);
   thrust_ /= PodParams::fric_thruts_to_thrust_ratio;
   thrust_ = std::min(thrust_, double(PodParams::max_thrust));
   thrust_ = std::max(thrust_, 0.0);
   thrust = round(thrust_);
   return thrust;
}

class Planner3: public Planner {
public:
   virtual void plan(icoord target, Runner *runner_) {
      runner = runner_;
      gt_regulator->reset(ArenaParams::station_tolerance, runner);
   }

   virtual instruction act(void);

private:
   unique_ptr<GoToTargetRegulator> gt_regulator = make_unique<GoToTargetRegulator>();
};

instruction Planner3::act(void)
{
   auto target = to_dcoord(runner->station(0));
   auto pos = to_dcoord(runner->cpos());
   //auto angle = angleabs2angle(runner->cangle(), target, pos);
   auto [ tc, thrust, boost, shield ] = gt_regulator->act(target);
   if (gt_regulator->is_pointing) {
      auto v1 = runner->vels.back().norm();
      if (v1 > 10.) {
         auto s = v1 / (1.0 - ArenaParams::friction_fac);
         auto num_turns = s / v1 * 0.5;
         auto dist = pos.dist_pnts(target);
         auto turns_to_targ = dist / v1;
         if (turns_to_targ <= num_turns)
            return make_tuple(runner->station(1), 0, false, false);
      }
   }
   return make_tuple(tc, thrust, boost, shield);
}

class Arena {
public:
    Arena(const string &name_, const list<icoord> &cps, const ArenaParams &a_params) :
        name(name_), stations(cps), params(a_params) {}

    bool pointInTrack(const icoord &p) const;

    string name;
    list <icoord> stations;
    ArenaParams params;
};

bool Arena::pointInTrack(const icoord &p) const
{
   for (auto &c: stations)
   {
      if (fabs(c.x - p.x) < 70 && fabs(c.y - p.y) < 70)
         return true;
   }
   return false;
}

class ArenaDetector{
public:
    ArenaDetector(void);

    const Arena *detect(const list<icoord> &stations);

    list<Arena> tracks;
    const Arena *detected_track = nullptr;
    list<icoord> stations;
};


ArenaDetector::ArenaDetector(void)
{
    tracks.push_back({"hostile",    {  { 13890, 1958 }, { 8009,  3263 },  { 2653,  7002 },                             { 10035, 5969           }                           },  HostileParams()          });
    tracks.push_back({"hostile2",   {  { 9409,  7247 }, { 5984,  4264 },  { 14637, 1420 },                             { 3470,  7203 }                           },            HostileParams()          });
    tracks.push_back({"pyramid",    {  { 7637,  5988 }, { 3133,  7544 },  { 9544,  4408 },                             { 14535, 7770 },                                        { 6312,  4294             },                 { 7782,  851  }                   }, PyramidParams()    });
    tracks.push_back({"triangle",   {  { 6012,  5376 }, { 11308, 2847 },  { 7482,  6917 }                           }, TriangleParams()         });
    tracks.push_back({"dalton",     {  { 9589,  1395 }, { 3618,  4434 },  { 8011,  7920 },                             { 13272, 5530            }                           }, DaltonParams()           });
    tracks.push_back({"makbilit",   {  { 12942, 7222 }, { 5655,  2587 },  { 4107,  7431 },                             { 13475, 2323 }                           },            MacbilitParams()         });
    tracks.push_back({"arrow",      {  { 10255, 4946 }, { 6114,  2172 },  { 3048,  5194 },                             { 6276,  7762 },                                        { 14119, 7768 },                             { 13897, 1216 }                   }, ArrowParams()      });
    tracks.push_back({"Shosh",      {  { 9116,  1857 }, { 5007,  5288 },  { 11505, 6074 }                           }, ShoshParams()            });
    tracks.push_back({"Til",        {  { 10558, 5973 }, { 3565,  5194 },  { 13578, 7574 },                             { 12430, 1357 }                           },            TilParams()              });
    tracks.push_back({"trapez",     {  { 11201, 5443 }, { 7257,  6657 },  { 5452,  2829 },                             { 10294, 3341 }                           },            TrapezParams()           });
    tracks.push_back({"Mehumash",   {  { 4049,  4630 }, { 13054, 1928 },  { 6582,  7823 },                             { 7494,  1330 },                                        { 12701, 7080 }                           }, MehumashParams() });
    tracks.push_back({"Trampoline", {  { 3307,  7251 }, { 14572, 7677 },  { 10588, 5072 },                             { 13100, 2343 },                                        { 4536,  2191 },                             { 7359,  4930 }                   }, TrampolineParams() });
    tracks.push_back({"Zigzag",     {  { 10660, 2295 }, { 8695,  7469 },  { 7225,  2174 },                             { 3596,  5288 },                                        { 13862, 5092 }                           }, ZigzagParams()   });
}

const Arena * ArenaDetector::detect(const list<icoord> &stations_)
{
   int num_tracks = 0;
   stations = stations_;
   for (const auto &track: tracks)
   {
      bool detected = true;
      for (const auto &cp: stations_)
         if (!track.pointInTrack(cp)) {
            detected = false;
            break;
         }
      if (detected) {
         ++num_tracks;
         detected_track = &track;
      }
   }
   if (num_tracks == 1) {
      cerr << "Single Arena detected: " << detected_track->name << endl;
      return detected_track;
   }
   return nullptr;
}


vector<unique_ptr<Runner> > Runner::me;
vector<unique_ptr<Runner> > Runner::other;

class Defender: public Runner
{
public:
   Defender(int num_laps, bool is_me): Runner(num_laps, is_me)
   {


   }
};



int main()
{
#ifndef TEST
    int num_laps;
    cin >> num_laps; cin.ignore();
    Runner::push_me(make_unique<Runner>(num_laps, true));
    Runner::push_me(make_unique<Defender>(num_laps, true));
    Runner::push_opponent(make_unique<Runner>(num_laps, false));
    Runner::push_opponent(make_unique<Runner>(num_laps, false));


    int check_point_count;
    cin >> check_point_count; cin.ignore();
    std::list<icoord> stations;
    for (int i = 0; i < check_point_count; i++) {
        int xs;
        int ys;
        cin >> xs >> ys; cin.ignore();
        stations.push_back({xs, ys});
    }
    ArenaDetector detector;
    if (!detector.detect(stations)) {
       cerr << " ******COULD NOT DETECT ARENA VERY BAD *********" << endl;
       exit(255);
    }


    // game loop
    while (true) {
        for (int i = 0; i < 2; i++) {
            int x; // x position of your pod
            int y; // y position of your pod
            int vx; // x speed of your pod
            int vy; // y speed of your pod
            int angle; // angle of your pod
            int nextCheckPointId; // next check point id of your pod
            cin >> x >> y >> vx >> vy >> angle >> nextCheckPointId; cin.ignore();
        }
        for (int i = 0; i < 2; i++) {
            int x2; // x position of the opponent's pod
            int y2; // y position of the opponent's pod
            int vx2; // x speed of the opponent's pod
            int vy2; // y speed of the opponent's pod
            int angle2; // angle of the opponent's pod
            int nextCheckPointId2; // next check point id of the opponent's pod
            cin >> x2 >> y2 >> vx2 >> vy2 >> angle2 >> nextCheckPointId2; cin.ignore();
        }

        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;


        // You have to output the target position
        // followed by the power (0 <= power <= 200)
        // i.e.: "x y power"
        cout << "8000 4500 100" << endl;
        cout << "8000 4500 100" << endl;
    }

#else
    test_line_intersect();
    test_rotate();
    test_circ_from_2points_and_rad();
    test_location_along_segment();
    test_angle_between();
    test_angleabs2angle();
#endif
    return 0;
}
