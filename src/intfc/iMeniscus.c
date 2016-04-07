/* This is where the meniscus will go */

/* Meniscus Helper Function*/

EXPORT double meniscus_profile_func(
        double *coords
        )
{
    int dim;
    dim = meniscus_params->dim;
    // It is a plane at z = z0 and at the boundary state, it's becoming a ramp
    double dist;
    // if x coordinate is less than certain number from left boundary state, dist = ramp
    // if x from right boundary state, dist = ramp
    // otherwise, dist = z0 - coords[dim-1]
}
