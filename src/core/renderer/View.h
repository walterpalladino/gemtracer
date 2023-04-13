#pragma once

#include "graphics/ImageBuffer.h"
#include "core/renderer/Camera.h"

class View
{
private:
	ImageBuffer imageBuffer;
	Camera camera;

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

	void RenderScene(float fAngle,
					 float fZPos,
					 float fXPos,
					 float fYPos);
	ImageBuffer GetImageBuffer(void) { return imageBuffer; }
};
