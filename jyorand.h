/*
		使用时先new一个对象：Rand jyorandengine;
		接着调用方法：jyorandengine.jyoRandGetInteger<int>(0, 100); 
*/

#ifndef JYORAND_H
#define JYORAND_H

#include <iostream>
#include <random>
#include <ctime>
#include <chrono>

class Rand
{
		private:
		std::mt19937 randEngineFunc[3];
		time_t seed[3] = {
			std::random_device {}(), 
			time(0), 
			std::chrono::system_clock::now().time_since_epoch().count()
		};
		int lasttype = 1;
		
		int typeChoose()
		{
			int type = 1;
			for(int i = 0; i <= 4; i ++)
			{
				static std::uniform_int_distribution<int> tp(0, 100);
				type = tp(randEngineFunc[type])%3;
				static std::bernoulli_distribution bd;
				if(bd(randEngineFunc[type])==true)
					i --;
			}
			if(type==lasttype)
			{
				type ++;
				type %= 3;
			}
			return lasttype = type;
		}
		public:
		//构造函数 
		Rand()
		{
			for(int i = 0;i <= 2;i ++)
				randEngineFunc[i].seed(seed[i]);
		}
		//均匀分布，[a,b]范围的整数 
		template <typename T>
		T jyoRandGetInteger(T a, T b)
		{
			int type = typeChoose();
			std::uniform_int_distribution<T> d(a, b);
			return d(randEngineFunc[type]);
		}
		//均匀分布，[a,b]范围的实数 
		template <typename T>
		T jyoRandGetReal(T a, T b)
		{
			int type = typeChoose();
			std::uniform_real_distribution<T> d(a, b);
			return d(randEngineFunc[type]);
		}
		//p概率返回true 
		bool jyoRandGetBool(double p = 0.5)
		{
			int type = typeChoose();
			std::bernoulli_distribution d(p);
			return d(randEngineFunc[type]);
		}
		//均值m，标准差s的正态分布结果 
		template <typename T>
		T jyoRandGetNormal(T m = 0.0, T s = 1.0)
		{
			int type = typeChoose();
			std::normal_distribution<T> d(m,s);
			return d(randEngineFunc[type]);
		}
};
#endif