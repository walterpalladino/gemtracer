
#include <math.h>

#include "math/MathCommons.h"
#include "math/Vector3d.h"
#include "core/materials/Bump.h"

void get_Bump_Normal(TWAVE *bump, Vector3d POINT, Vector3d *R)
{

    Vector3d diff;
    double dist;
    double amp;

    //    subV( POINT, bump->center, &diff );
    diff = POINT - bump->center;
    //    normalize( diff, &diff );
    diff.Normalize();
    //    dist = lenV( POINT, bump->center );
    dist = (POINT - bump->center).Magnitude();

    dist /= bump->wavelength;
    dist += bump->phase;

    if (bump->attenuation < 1.0)
        amp = bump->amplitude * pow(bump->attenuation, dist);
    else
        amp = bump->amplitude;

    amp = cos(dist * PI * 2.0);

    //	mulS ( diff, amp, R ) ;
    *R = diff * amp;

    //    normalize( *R, R );

    /*
        subV( POINT, bump->center, &diff );
        dist = lenV( POINT, bump->center );

        dist /= bump->wavelength;
        dist += bump->phase;

        if (dist == 0.0)
          dist = 1.0;

        amp = bump->amplitude * sin ( dist * PI * 4.0 ) ;

        mulS ( diff, amp/bump->wavelength, R ) ;
        normalize( *R, R );
    */
}
