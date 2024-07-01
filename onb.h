#ifndef ONBH
#define ONBH

#include <iostream>
#include "vec3.h"

class onb{
    public:
        vec3 axis[3];
    
    public:
    onb() {}
    // 以给定的一个向量为轴
    onb(const vec3& w){
        // 单位化给定的向量
        vec3 unit_w = unit_vector(w);
        // 找一个坐标轴，但不能和w平行
        vec3 a = (fabs(unit_w.x()) > 0.9 ? vec3(0, 1, 0) : vec3(1, 0, 0));
        // 叉乘得到一个垂直于这两个向量的坐标系的一轴
        vec3 v = unit_vector(cross(unit_w, a));
        // 再次叉乘得到坐标系的另一轴
        vec3 u = cross(unit_w, v);
        axis[0] = u, axis[1] = v, axis[2] = unit_w;
    }

    vec3 operator[](int i) const {return axis[i];}
    vec3& operator[](int i) {return axis[i];}

    vec3 u() const {return axis[0];}
    vec3 v() const {return axis[1];}
    vec3 w() const {return axis[2];}

    // 从全局坐标转换为局部坐标
    vec3 local(double a, double b, double c) const {
        return a*u() + b*v() + c*w();
    }

    vec3 local(const vec3& a) const {
        return a.x()*u() + a.y()*v() + a.z()*w();
    }
};

#endif