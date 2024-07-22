#ifndef PERLINH
#define PERLINH

#include "vec3.h"

// 前向声明
extern Rand jyorandengine;

// 线性插值
inline double trilinear_interp(vec3 c[2][2][2], double u, double v,
                               double w) {
  // hermite cubic舍入插值
  double uu = u * u * (3 - 2 * u);
  double vv = v * v * (3 - 2 * v);
  double ww = w * w * (3 - 2 * w);

  double accum = 0;
  for (int i = 0; i <= 1; i++)
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++) {
        vec3 weight_v(u-i, v-j, w-k);
        accum += (i * uu + (1 - i) * (1 - uu)) * (j * vv + (1 - j) * (1 - vv)) *
                 (k * ww + (1 - k) * (1 - ww)) * dot(c[i][j][k], weight_v);
      }

  return accum;
}

class perlin {
 private:
  std::vector<vec3> ranvec3;
  std::vector<int> permx;
  std::vector<int> permy;
  std::vector<int> permz;
  std::vector<vec3> perlin_generate();
  void permute(std::vector<int>& v);
  std::vector<int> perlin_generate_perm();

 public:
  perlin() {
    ranvec3 = perlin_generate();
    permx = perlin_generate_perm();
    permy = perlin_generate_perm();
    permz = perlin_generate_perm();
  }
  // 输入一个点坐标，返回一个[0, 1]的随机浮点数值
  double noise(const vec3& p) const {
    // floor函数向下取整，uvw都是[0,1)的浮点数
    double u = p.x() - floor(p.x());
    double v = p.y() - floor(p.y());
    double w = p.z() - floor(p.z());

    // 向下取整的该点近似坐标
    int i = floor(p.x());
    int j = floor(p.y());
    int k = floor(p.z());

    vec3 c[2][2][2];
    for (int di = 0; di <= 1; di++)
      for (int dj = 0; dj <= 1; dj++)
        for (int dk = 0; dk <= 1; dk++)
          // 异或运算合并哈希
          c[di][dj][dk] =
              ranvec3[permx[(i + di) & 255] ^ permy[(j + dj) & 255] ^
                        permz[(k + dk) & 255]];
    return trilinear_interp(c, u, v, w);
  }
};

// 随机获得一个长度为256的，装满了[0,1]浮点数的浮点向量数组
std::vector<vec3> perlin::perlin_generate() {
  std::vector<vec3> v(256);
  for (int i = 0; i < 256; i++)
    v[i] =
        unit_vector(vec3(-1 + 2 * jyorandengine.jyoRandGetReal<double>(0, 1),
                         -1 + 2 * jyorandengine.jyoRandGetReal<double>(0, 1),
                         -1 + 2 * jyorandengine.jyoRandGetReal<double>(0, 1)));
  return v;
}

// 对每一个元素进行数组内元素随机交换
void perlin::permute(std::vector<int>& v) {
  for (int i = v.size() - 1; i > 0; i--) {
    int target = int(jyorandengine.jyoRandGetReal<double>(0, 1) * (double)(i + 1));
    std::swap(v[i], v[target]);
  }
}

// 获得一个长度为256的排列
std::vector<int> perlin::perlin_generate_perm() {
  std::vector<int> v(256, 0);
  for (int i = 0; i < 256; i++) v[i] = i;
  permute(v);
  return v;
}

#endif