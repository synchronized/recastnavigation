//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#ifndef DEBUGDRAW_H
#define DEBUGDRAW_H

// Some math headers don't have PI defined.
static const float DU_PI = 3.14159265f;

enum duDebugDrawPrimitives
{
	DU_DRAW_POINTS,
	DU_DRAW_LINES,
	DU_DRAW_TRIS,
	DU_DRAW_QUADS
};

/// Abstract debug draw interface.
/// 抽象调试绘制接口
struct duDebugDraw
{
	virtual ~duDebugDraw() = 0;

	virtual void depthMask(bool state) = 0;

	virtual void texture(bool state) = 0;

	/// Begin drawing primitives.
	///  @param prim [in] primitive type to draw, one of rcDebugDrawPrimitives.
	///  @param size [in] size of a primitive, applies to point size and line width only.
	/// 开始绘制基元
	///  @param prim [in] 绘制的基元类型, 必须是枚举duDebugDrawPrimitives中的一个
	///  @param size [in] 图元的大小，仅适用于点大小和线宽
	virtual void begin(duDebugDrawPrimitives prim, float size = 1.0f) = 0;

	/// Submit a vertex
	///  @param pos [in] position of the verts.
	///  @param color [in] color of the verts.
	/// 提交一个顶点
	///  @param pos [in] 顶点位置
	///  @param color [in] 顶点颜色
	virtual void vertex(const float* pos, unsigned int color) = 0;

	/// Submit a vertex
	///  @param x,y,z [in] position of the verts.
	///  @param color [in] color of the verts.
	/// 提交一个顶点
	///  @param x,y,z [in] 顶点位置
	///  @param color [in] 顶点颜色
	virtual void vertex(const float x, const float y, const float z, unsigned int color) = 0;

	/// Submit a vertex
	///  @param pos [in] position of the verts.
	///  @param color [in] color of the verts.
	///  @param uv [in] the uv coordinates of the verts.
	/// 提交一个顶点
	///  @param pos [in] 顶点位置
	///  @param color [in] 顶点颜色
	///  @param uv [in] 顶点的uv坐标
	virtual void vertex(const float* pos, unsigned int color, const float* uv) = 0;

	/// Submit a vertex
	///  @param x,y,z [in] position of the verts.
	///  @param color [in] color of the verts.
	///  @param u,v [in] the uv coordinates of the verts.
	/// 提交一个顶点
	///  @param x,y,z [in] 顶点位置
	///  @param color [in] 顶点颜色
	///  @param u,v [in] 顶点的uv坐标
	virtual void vertex(const float x, const float y, const float z, unsigned int color, const float u, const float v) = 0;

	/// End drawing primitives.
	/// 结束绘制基元
	virtual void end() = 0;

	/// Compute a color for given area.
	/// 在区域中计算颜色
	virtual unsigned int areaToCol(unsigned int area);
};

inline unsigned int duRGBA(int r, int g, int b, int a)
{
	return ((unsigned int)r) | ((unsigned int)g << 8) | ((unsigned int)b << 16) | ((unsigned int)a << 24);
}

inline unsigned int duRGBAf(float fr, float fg, float fb, float fa)
{
	unsigned char r = (unsigned char)(fr*255.0f);
	unsigned char g = (unsigned char)(fg*255.0f);
	unsigned char b = (unsigned char)(fb*255.0f);
	unsigned char a = (unsigned char)(fa*255.0f);
	return duRGBA(r,g,b,a);
}

unsigned int duIntToCol(int i, int a);
void duIntToCol(int i, float* col);

//颜色乘法
inline unsigned int duMultCol(const unsigned int col, const unsigned int d)
{
	const unsigned int r = col & 0xff;
	const unsigned int g = (col >> 8) & 0xff;
	const unsigned int b = (col >> 16) & 0xff;
	const unsigned int a = (col >> 24) & 0xff;
	return duRGBA((r*d) >> 8, (g*d) >> 8, (b*d) >> 8, a);
}

//颜色变暗
inline unsigned int duDarkenCol(unsigned int col)
{
	return ((col >> 1) & 0x007f7f7f) | (col & 0xff000000);
}

//颜色从ca到cb渐变的插值,u代表颜色变化的程度
inline unsigned int duLerpCol(unsigned int ca, unsigned int cb, unsigned int u)
{
	const unsigned int ra = ca & 0xff;
	const unsigned int ga = (ca >> 8) & 0xff;
	const unsigned int ba = (ca >> 16) & 0xff;
	const unsigned int aa = (ca >> 24) & 0xff;
	const unsigned int rb = cb & 0xff;
	const unsigned int gb = (cb >> 8) & 0xff;
	const unsigned int bb = (cb >> 16) & 0xff;
	const unsigned int ab = (cb >> 24) & 0xff;

	unsigned int r = (ra*(255-u) + rb*u)/255;
	unsigned int g = (ga*(255-u) + gb*u)/255;
	unsigned int b = (ba*(255-u) + bb*u)/255;
	unsigned int a = (aa*(255-u) + ab*u)/255;
	return duRGBA(r,g,b,a);
}

//设置颜色的透明度
inline unsigned int duTransCol(unsigned int c, unsigned int a)
{
	return (a<<24) | (c & 0x00ffffff);
}


// 计算长方体颜色
void duCalcBoxColors(unsigned int* colors, unsigned int colTop, unsigned int colSide);

// 绘制圆柱体
void duDebugDrawCylinderWire(struct duDebugDraw* dd, float minx, float miny, float minz,
							 float maxx, float maxy, float maxz, unsigned int col, const float lineWidth);

// 绘制正方体(12条棱)
void duDebugDrawBoxWire(struct duDebugDraw* dd, float minx, float miny, float minz,
						float maxx, float maxy, float maxz, unsigned int col, const float lineWidth);

// 绘制弧度
void duDebugDrawArc(struct duDebugDraw* dd, const float x0, const float y0, const float z0,
					const float x1, const float y1, const float z1, const float h,
					const float as0, const float as1, unsigned int col, const float lineWidth);

// 绘制一个箭头
void duDebugDrawArrow(struct duDebugDraw* dd, const float x0, const float y0, const float z0,
					  const float x1, const float y1, const float z1,
					  const float as0, const float as1, unsigned int col, const float lineWidth);

// 绘制一个圆(40条边)
void duDebugDrawCircle(struct duDebugDraw* dd, const float x, const float y, const float z,
					   const float r, unsigned int col, const float lineWidth);

// 绘制一个叉
void duDebugDrawCross(struct duDebugDraw* dd, const float x, const float y, const float z,
					  const float size, unsigned int col, const float lineWidth);

// 绘制长方体(6个面)
void duDebugDrawBox(struct duDebugDraw* dd, float minx, float miny, float minz,
					float maxx, float maxy, float maxz, const unsigned int* fcol);

// 绘制圆柱体(三角形)
void duDebugDrawCylinder(struct duDebugDraw* dd, float minx, float miny, float minz,
						 float maxx, float maxy, float maxz, unsigned int col);

// 绘制一个xz坐标的格子
//  @param ox,oy,oz 起始点
//  @param w 列的数量(z轴)
//  @param h 行的数量(x轴)
//  @param size 格子的边长
//  @param col 颜色
//  @param lineWidth 线的宽度
void duDebugDrawGridXZ(struct duDebugDraw* dd, const float ox, const float oy, const float oz,
					   const int w, const int h, const float size,
					   const unsigned int col, const float lineWidth);


// Versions without begin/end, can be used to draw multiple primitives.
// 没有begin/end的版本, 能够用来进行多次绘制基元

// 绘制圆柱体(16条棱/显示4条棱)
void duAppendCylinderWire(struct duDebugDraw* dd, float minx, float miny, float minz,
						  float maxx, float maxy, float maxz, unsigned int col);

// 绘制长方体(线)
void duAppendBoxWire(struct duDebugDraw* dd, float minx, float miny, float minz,
					 float maxx, float maxy, float maxz, unsigned int col);

// 绘制长方体(点)
void duAppendBoxPoints(struct duDebugDraw* dd, float minx, float miny, float minz,
					   float maxx, float maxy, float maxz, unsigned int col);

// 绘制弧度
void duAppendArc(struct duDebugDraw* dd, const float x0, const float y0, const float z0,
				 const float x1, const float y1, const float z1, const float h,
				 const float as0, const float as1, unsigned int col);

// 绘制箭头
void duAppendArrow(struct duDebugDraw* dd, const float x0, const float y0, const float z0,
				   const float x1, const float y1, const float z1,
				   const float as0, const float as1, unsigned int col);

// 绘制圆形
void duAppendCircle(struct duDebugDraw* dd, const float x, const float y, const float z,
					const float r, unsigned int col);

// 绘制一个叉
void duAppendCross(struct duDebugDraw* dd, const float x, const float y, const float z,
				   const float size, unsigned int col);

// 绘制长方体(6个面)
void duAppendBox(struct duDebugDraw* dd, float minx, float miny, float minz,
				 float maxx, float maxy, float maxz, const unsigned int* fcol);

// 绘制圆柱体(三角形)
void duAppendCylinder(struct duDebugDraw* dd, float minx, float miny, float minz,
					  float maxx, float maxy, float maxz, unsigned int col);


class duDisplayList : public duDebugDraw
{
	float* m_pos;
	unsigned int* m_color;
	int m_size;
	int m_cap;

	duDebugDrawPrimitives m_prim;
	float m_primSize;
	bool m_depthMask;

	void resize(int cap);

public:
	duDisplayList(int cap = 512);
	virtual ~duDisplayList();
	virtual void depthMask(bool state);
	virtual void begin(duDebugDrawPrimitives prim, float size = 1.0f);
	virtual void vertex(const float x, const float y, const float z, unsigned int color);
	virtual void vertex(const float* pos, unsigned int color);
	virtual void end();
	void clear();
	void draw(struct duDebugDraw* dd);
private:
	// Explicitly disabled copy constructor and copy assignment operator.
	duDisplayList(const duDisplayList&);
	duDisplayList& operator=(const duDisplayList&);
};


#endif // DEBUGDRAW_H
