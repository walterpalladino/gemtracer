#pragma once

#include "graphics/ImageBuffer.h"
#include "core/renderer/Camera.h"
#include "math/Vector3d.h"

class View
{
private:
	ImageBuffer imageBuffer;
	Camera camera;

	Vector3d TraceRay(float x, float y);

public:
	long width;
	long height;

	View(void);

	static View *GetInstance()
	{
		static View instance;
		return &instance;
	};

	long Init(void);
	long Init(int width, int height);

	void SetCamera(Vector3d location,
				   Vector3d direction,
				   Vector3d right,
				   Vector3d up,
				   double angle);

	void LoadScene(char *sceneName);

	void RenderScene(bool useSSAA);
	ImageBuffer GetImageBuffer(void) { return imageBuffer; }
};
