// Copyright (c) 2016, Michael P. Howard. All rights reserved.
// This software is released under the BSD 3-Clause license.

#ifndef SOPHIA_VECTORMATH_H_
#define SOPHIA_VECTORMATH_H_


//! Two-element vector
template <typename T>
struct vec2
    {
    T x,y;
    vec2() : x(0), y(0) {}
    vec2(T _x, T _y) : x(_x), y(_y) {}
    
    vec2<T> operator+ (const vec2<T>& a) const
        {
	    return vec2<T>(x+a.x,y+a.y);
        }
	
    vec2<T> operator- (const vec2<T>& a) const
        {
	    return vec2<T>(x-a.x,y-a.y);
        }
	
	vec2<T>& operator+= (const vec2<T>& a)
        {
	    x += a.x;
	    y += a.y;
	    return *this;
        }
        
	vec2<T>& operator-= (const vec2<T>& a)
        {
	    x -= a.x;
	    y -= a.y;
	    return *this;
        }
    
    bool operator== (const vec2<T>& a) const
        {
        return (x == a.x && y == a.y);
        }
        
    bool operator!= (const vec2<T>& a) const
        {
        return (x != a.x || y != a.y);
        }
    
    double dot(const vec2<T>& a) const
        {
        return x*a.x + y*a.y;
        }
        
    };
typedef vec2<float> float2;
typedef vec2<double> double2;

//! Three-element vector
template <typename T>
struct vec3
    {
    T x,y,z;
    vec3() : x(0), y(0), z(0) {}
    vec3(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
    
    vec3<T> operator+ (const vec3<T>& a) const
        {
	    return vec3<T>(x+a.x,y+a.y,z+a.z);
        }
	
    vec3<T> operator- (const vec3<T>& a) const
        {
	    return vec3<T>(x-a.x,y-a.y,z-a.z);
        }
	
	vec3<T>& operator+= (const vec3<T>& a)
        {
	    x += a.x;
	    y += a.y;
	    z += a.z;
	    return *this;
        }
        
	vec3<T>& operator-= (const vec3<T>& a)
        {
	    x -= a.x;
	    y -= a.y;
	    z -= a.z;
	    return *this;
        }
    
    bool operator== (const vec3<T>& a) const
        {
        return (x == a.x && y == a.y && z == a.z);
        }
        
    bool operator!= (const vec3<T>& a) const
        {
        return (x != a.x || y != a.y || z != a.z);
        }
    
    double dot(const vec3<T>& a) const
        {
        return x*a.x + y*a.y + z*a.z;
        }
        
    //! this.cross(a) = this x a
    vec3<T> cross(const vec3<T>& a) const
        {
        return vec3<T>(y*a.z-z*a.y, z*a.x-x*a.z, x*a.y-y*a.x);
        }
    };
typedef vec3<float> float3;
typedef vec3<double> double3;

//! Four-element vector
template <typename T>
struct vec4
    {
    T x,y,z,n;
    vec4() : x(0), y(0), z(0), n(0) {}
    vec4(T _x, T _y, T _z, T _n) : x(_x), y(_y), z(_z), n(_n){}
    
    vec4<T> operator+ (const vec4<T>& a) const
        {
	    return vec4<T>(x+a.x,y+a.y,z+a.z,n+a.n);
        }
	
    vec4<T> operator- (const vec4<T>& a) const
        {
	    return vec3<T>(x-a.x,y-a.y,z-a.z,n-a.n);
        }
	
	vec4<T>& operator+= (const vec4<T>& a)
        {
	    x += a.x;
	    y += a.y;
	    z += a.z;
		n += a.n;
	    return *this;
        }
        
	vec4<T>& operator-= (const vec4<T>& a)
        {
	    x -= a.x;
	    y -= a.y;
	    z -= a.z;
		n -= a.n;
	    return *this;
        }
    
    bool operator== (const vec4<T>& a) const
        {
        return (x == a.x && y == a.y && z == a.z && n == a.n);
        }
        
    bool operator!= (const vec4<T>& a) const
        {
        return (x != a.x || y != a.y || z != a.z || z != a.n);
        }
    
    double dot(const vec4<T>& a) const
        {
        return x*a.x + y*a.y + z*a.z + n*a.n;
        }
        
    };
typedef vec4<float> float4;
typedef vec4<double> double4;
#endif // SOPHIA_VECTORMATH_H_
