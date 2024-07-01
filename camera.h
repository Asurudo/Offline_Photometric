#ifndef CAEERAH
#define CAMERAH

#include "jyorand.h"
#include "ray.h"

extern vec3 randomInUnitDisk();
extern Rand jyorandengine;

class camera {
 public:
  // 相机所在位置
  vec3 origin;
  // 画布左下角坐标
  vec3 lower_left_corner;
  // 画布宽度
  vec3 horizontal;
  // 画布高度
  vec3 vertical;
  // uvw是以相机为原点的标准向量
  vec3 u, v, w;
  // 镜头半径
  double lensRadius;
  // 快门的开关时间
  double time0, time1;

  // lookfrom是相机位置，也就是origin
  // lookat是相机的朝向位置
  // vfov是相机在垂直方向上从屏幕顶端扫描到屏幕底端岔开的视角角度，设为90°，也就是二分之PI
  // aspect是屏幕的宽高比，设为2:1
  // aperture是光圈的直径
  // focus是焦距，也就是光圈到成像平面的距离
  camera(vec3 lookfrom, vec3 lookat, double vfov, double aspect, double t0,
         double t1, double aperture = 0.0, double focus = 1.0) {
    time0 = t0;
    time1 = t1;

    // vup是一个固定值,是一个指向世界正上方的单位向量，用于辅助运算
    const vec3 vup = vec3(0, 1, 0);
    // theta为视角范围
    double theta = vfov * PI / 180;
    // 画布的半高
    double halfHeight = tan(theta / 2.0) * focus;
    // 画布的半宽
    double halfWidth = aspect * halfHeight;

    origin = lookfrom;
    lensRadius = aperture / 2.0;

    w = unit_vector(lookfrom - lookat);
    u = unit_vector(cross(vup, w));
    v = cross(w, u);

    lower_left_corner = origin - halfWidth * u - halfHeight * v - focus * w;
    horizontal = 2 * halfWidth * u;
    vertical = 2 * halfHeight * v;
  }

  ray get_ray(double s, double t) {
    // 以（0，0，0）为原点的随机出来的点
    vec3 randomPoint = lensRadius * randomInUnitDisk();
    // 以相机坐标为原点随机出来的点
    vec3 offset = randomPoint.x() * u + randomPoint.y() * v;

    double randNum = jyorandengine.jyoRandGetNormal<double>((time1-time0)/2.0, (time1-time0)/5);
    if (randNum < time0)
      randNum = time0;
    else if (randNum > time1)
      randNum = time1;
    double time = randNum;
    // std::cout << randNum << std::endl;
    // 一条从相机到视点的射线，也就是一条光线
    return ray(
        origin + offset,
        lower_left_corner + s * horizontal + t * vertical - origin - offset,
        time);
  }
};

#endif