#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <memory>
#include <utility>

using namespace std;

#define PI 3.14159265

float radians(float ang) {
   return ang * PI / 180.;
}

float degrees(float rad) {
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

   Coord operator - (const Coord &c1) {
      return {x - c1.x, y - c1.y};
   }

   float norm(void) const {
      return sqrt(x*x + y*y);
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

   Coord operator / (T t) {
      return {x / t, y / t};
   }

   Coord operator * (T t) {
      return {x * t, y * t};
   }

   T cross(const Coord &c) {
      return x * c.y - y * c.x;
   }

   Coord unit_vec(void) {
      return (*this)/norm();
   }

   float angle_between(const Coord &v2) {
      auto v1_u = unit_vec();
      auto v2_u = v2.unit_vector();
      auto c = clip(dot(*this, v2), -1, 1);
      auto ret = degrees(arccos(c));
      auto cr = v1_u.cross(v2_u);
      return (cr >= 0) ? ret : -ret;
   }
};
template <typename T=int>
float relAngle(const Coord<T> &p1, const Coord<T> &p2) {
   return atan2(Coord<T>::det(p1, p2), Coord<T>::dot(p1, p2));
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

typedef Coord<int> icoord;
typedef Coord<float> FCoord;

typedef icoord ivec;

struct PodParams {};
struct ArenaParams {

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
   return ArenaParams();
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

class Runner;

class Predictor
{
public:
    Predictor(Runner *runner);

private:


};

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


/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/

class Runner
{
public:
   Runner(int num_laps, bool is_me) {


   }

   static void push_me(unique_ptr<Runner> r) {
       me.push_back(std::move(r));
   }

   static void push_opponent(unique_ptr<Runner> r) {
       other.push_back(std::move(r));
   }

private:
   static                  vector<unique_ptr<Runner> > me;
   static                  vector<unique_ptr<Runner> > other;

   int                     num_laps = 0;
   ArenaParams             arena_params;
   PodParams               pod_params;
   list<icoord>             pos;
   list<int>               angles;
   list<ivec>              vels;
   list<int>               thrusts;
   int                     next_cp = 0;
   int                     prev_cp = 0;
   list<icoord>             stations;
   int                     time_left_to_punch = 0;
   int                     passed_stations = -1;
   bool                    is_me = false;
   bool                    boost = false;
   int                     boost_turns = 0;
   shared_ptr<Predictor>   predictor;
   list<icoord>             predictions;
};

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

    return 0;
}
