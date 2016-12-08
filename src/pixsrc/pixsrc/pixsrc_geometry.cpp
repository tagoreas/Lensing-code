#ifndef PS_SIT
#ifdef INTEGER_TYPE_LL
typedef long long PS_SIT;
typedef unsigned long long PS_unsignedSIT;
#endif
#ifdef INTEGER_TYPE_L
typedef long PS_SIT;
typedef unsigned long PS_unsignedSIT;
#endif
#ifdef INTEGER_TYPE_I
typedef int PS_SIT;
typedef unsigned int PS_unsignedSIT;
#endif
#endif

#ifndef PS_FPT
#ifdef SINGLE_PRECISION
typedef float PS_FPT;
#endif
#ifdef DOUBLE_PRECISION
typedef double PS_FPT;
#endif
#endif



#include "pixsrc_memory_templates.cpp"
#include "pixsrc_geometry.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_constants.hpp"
#include <cmath>
#include <cstring>

void pixsrc_geometry::getcentroid(PS_SIT numverts, double *hull, double *ctrx, double *ctry)
{
    // have to compute area here b/c we need the signed area
    // function below gives unsigned area
    double area = 0;
    *ctrx = *ctry = 0;

    PS_SIT j = numverts - 1;
    for( PS_SIT i=0; i<numverts; ++i )
    {
        area+=(hull[j*2]+hull[i*2])*(hull[j*2+1]-hull[i*2+1]);
        *ctrx += ( hull[i*2  ]+hull[j*2  ] ) *
            ( hull[i*2]*hull[j*2+1] - hull[j*2]*hull[i*2+1] );
        *ctry += ( hull[i*2+1]+hull[j*2+1] ) *
            ( hull[i*2]*hull[j*2+1] - hull[j*2]*hull[i*2+1] );

        j=i;
    }

    area /= 2.0;
    *ctrx /= ( 6.0*area );
    *ctry /= ( 6.0*area );
}
double pixsrc_geometry::areatriangle(double *coord)
{
    double area=0;
    PS_SIT j=2;

    for(PS_SIT i=0; i<3; i++)
    {
        area+=(coord[j*2]+coord[i*2])*(coord[j*2+1]-coord[i*2+1]);
        j=i;
    }

    return std::fabs(area*0.5);
}

double pixsrc_geometry::areapoly(double **coord, PS_SIT numverts)
{
    if (numverts<3)
        return 0;

    double area=0;
    PS_SIT j=numverts-1;

    for(PS_SIT i=0; i<numverts; i++)
    {
        area+=(coord[j][0]+coord[i][0])*(coord[j][1]-coord[i][1]);
        j=i;
    }

    return std::fabs(area*0.5);
}
double pixsrc_geometry::areapoly(double *coord, PS_SIT numverts)
{
    if (numverts<3)
        return 0;

    double area=0;
    PS_SIT j=numverts-1;

    for(PS_SIT i=0; i<numverts; i++)
    {
        area+=(coord[2*j]+coord[2*i])*(coord[2*j+1]-coord[2*i+1]);
        j=i;
    }

    return std::fabs(area*0.5);
}
double pixsrc_geometry::getmaxdist(double *coord, PS_SIT numverts)
{
    double dwinning = 0;
    for(PS_SIT i=0; i<numverts; i++)
    {
        for(PS_SIT j=i+1; j<numverts; j++)
        {
            double d = OPERA distance(coord[i*2],coord[i*2+1],coord[j*2],coord[j*2+1]);
            if(d>dwinning)
                dwinning=d;
        }
    }
    return dwinning;
}
bool pixsrc_geometry::isintricross(double *pos)
{
    double cross[3];
    PS_SIT i = 2;
    for( PS_SIT j=0; j<3; ++j )
    {
        // calculating cross product of vector of position
        // vectors of one vertex and next
        cross[j] = pos[i*2]*pos[j*2+1] - pos[i*2+1]*pos[j*2];
        i = j;
    }

    // if cross products circulate in same direction or
    // x4,y4 is on a line segment
    return ( (cross[0]<0 && cross[1]<0 && cross[2]<0) ||
             (cross[0]>0 && cross[1]>0 && cross[2]>0) ||
             !cross[0] || !cross[1] || !cross[2]    )  ? 1 : 0;
}
bool pixsrc_geometry::isintri(double *tri, double x4, double y4)
{
    // point x4,y4 is the one being tested
    // shold use this functino over isinpoly if poly=triangle

    double pos[6];
    for( PS_SIT j=0; j<3; ++j )
    {
        // shifting to coord. sys. centered on x4,y4
        pos[j*2  ] = tri[j*2  ] - x4;
        pos[j*2+1] = tri[j*2+1] - y4;
    }

    return isintricross( pos );
}
bool pixsrc_geometry::isintripadded(double *tri0, double x4, double y4)
{
    // point x4,y4 is the one being tested
    // shold use this functino over isinpoly if poly=triangle

    double tri[2], pos[6], r;
    // centroid
    double ctr[2] = { (tri0[0]+tri0[2]+tri0[4])/3.0,(tri0[1]+tri0[3]+tri0[5])/3.0 };
    for( PS_SIT j=0; j<3; ++j )
    {
        // here i shift the triangle coordinates by a small amount
        // in a direction pointing from centroid to vertex

        // shifting to coord. sys. w. centroid at origin
        tri[0] = tri0[j*2  ] - ctr[0];
        tri[1] = tri0[j*2+1] - ctr[1];
        r = std::sqrt( tri[0]*tri[0] + tri[1]*tri[1] );
        // shifting tri. cooridnates and then shifting to coord. sys.
        // centered on x4,y4
        pos[j*2  ] = tri0[j*2  ] + tri[0]/r*CONSTANT smallnumber - x4;
        pos[j*2+1] = tri0[j*2+1] + tri[1]/r*CONSTANT smallnumber - y4;
    }

    return isintricross( pos );
}
bool pixsrc_geometry::isintri(double pos[2][3], double x4, double y4)
{
    // point x4,y4 is the one being tested
    // shold use this functino over isinpoly if poly=triangle

    double pos2[6];
    for( PS_SIT j=0; j<3; ++j )
    {
        // shifting to coord. sys. centered on x4,y4
        pos2[j*2  ] = pos[0][j] - x4;
        pos2[j*2+1] = pos[1][j] - y4;
    }

    return isintricross( pos2 );
}
bool pixsrc_geometry::isintripadded(double pos[2][3], double x4, double y4)
{
    // point x4,y4 is the one being tested
    // shold use this functino over isinpoly if poly=triangle

    double tri[2], pos2[6], r;
    // centroid
    double ctr[2] = { (pos[0][0]+pos[0][1]+pos[0][2])/3.0,(pos[1][0]+pos[1][1]+pos[1][2])/3.0 };
    for( PS_SIT j=0; j<3; ++j )
    {
        // here i shift the triangle coordinates by a small amount
        // in a direction pointing from centroid to vertex

        // shifting to coord. sys. w. centroid at origin
        tri[0] = pos[0][j] - ctr[0];
        tri[1] = pos[1][j] - ctr[1];
        r = std::sqrt( tri[0]*tri[0] + tri[1]*tri[1] );
        // shifting tri. cooridnates and then shifting to coord. sys.
        // centered on x4,y4
        pos2[j*2  ] = pos[0][j] + tri[0]/r*CONSTANT smallnumber - x4;
        pos2[j*2+1] = pos[1][j] + tri[1]/r*CONSTANT smallnumber - y4;
    }

    return GEOM isintricross( pos2 );
}
bool pixsrc_geometry::isinpolycalc( double *p2, PS_SIT numv )
{
    PS_SIT y = numv-1;
    bool returner = 0;
    for(PS_SIT x=0; x<numv; x++)
    {
        // i forget how this works exactly
        if( ( (p2[x*2+1]<0 && p2[y*2+1]>=0) || (p2[y*2+1]<0 && p2[x*2+1]>=0) ) &&
            (p2[x*2]-p2[x*2+1]/(p2[y*2+1]-p2[x*2+1])*(p2[y*2]-p2[x*2]) < 0) )
            returner=!returner;
        y=x;
    }

    return returner;
}
bool pixsrc_geometry::isinpoly(double x0, double y0, double *p, PS_SIT numpoly)
{
    double p2[numpoly*2];
    for(PS_SIT x=0; x<numpoly; x++)
    {
        // shifting to coord. sys. centered on x0,y0
        p2[x*2  ] = p[x*2  ] - x0;
        p2[x*2+1] = p[x*2+1] - y0;
    }

    return isinpolycalc( p2,numpoly );
}
bool pixsrc_geometry::isinpolypadded(double x0, double y0, double *p, PS_SIT numpoly)
{
    double p2[numpoly*2], ctr[2], pos[2], r;
    GEOM getcentroid( numpoly,p,ctr,ctr+1 );
    for(PS_SIT x=0; x<numpoly; x++)
    {
        // here i shift the polygon coordinates by a small amount
        // in a direction pointing from centroid to vertex

        // shifting to coord. sys. w. centroid at origin
        pos[0] = p[x*2  ] - ctr[0];
        pos[1] = p[x*2+1] - ctr[1];
        r = std::sqrt( pos[0]*pos[0] + pos[1]*pos[1] );
        // shifting polygon cooridnates and then shifting to coord. sys.
        // centered on x0,y0
        p2[x*2  ] = p[x*2  ] + pos[0]/r*CONSTANT smallnumber - x0;
        p2[x*2+1] = p[x*2+1] + pos[1]/r*CONSTANT smallnumber - y0;
    }

    return isinpolycalc( p2,numpoly );
}
bool pixsrc_geometry::isonlinesegment(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double cross = (x1-x3)*(y2-y3)-(y1-y3)*(x2-x3);
    return ( OPERA equalszero(cross) && (x3-x1)*(x3-x2)<=0 && (y3-y1)*(y3-y2)<=0 ) ? 1 : 0;
}

double pixsrc_geometry::angle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    // gets angle between two vectors, starting and ending at ends of vectors. points 1 and 2 specify ends of vectors.
    // meaning .. 1 and 3 are ends and 2 is vertex
    // returns something between 0 180

    return acos( (x1*x3+y1*y3-x1*x2-y1*y2-x2*x3-y2*y3+x2*x2+y2*y2) /
                 sqrt( ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) *
                       ((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)) ) );
}
bool pixsrc_geometry::allindiffquadrants(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    double array[8] = {x1,y1,x2,y2,x3,y3,x4,y4};

    for(PS_SIT u=0; u<4; u++)
    {
        if( OPERA equalszero( array[2*u]) )
        {
            // if point is at center
            if(array[2*u+1]==0)
                return 0;
            array[2*u] = 0;

            if(array[2*u+1]>0)
                array[2*u] =  1;
            else
                array[2*u] = -1;
        }
        // don't need to test if previous has already been tested
        else if( OPERA equalszero(array[2*u+1]) )
        {
            array[2*u+1] = 0;

            if(array[2*u]>0)
                array[2*u+1] = -1;
            else
                array[2*u+1] =  1;
        }
    }

    for(PS_SIT x = 0; x < 4; x++)
        for(PS_SIT y = x+1; y < 4; y++)
            if(array[2*x]*array[2*y]>0 && array[2*x+1]*array[2*y+1]>0)
                return 0;
    return 1;
}
void pixsrc_geometry::norm_parallelogram_ic(double a, double b, double c, double d, double e, double f, double g, double h, bool isintriangle0, bool isintriangle1, bool isintriangle2, bool isintriangle3, PS_SIT *result)
{
    //if the "top" and "bottom" angles sum to larger than 180
    if( (angle(c,d,a,b,g,h)+angle(c,d,e,f,g,h)>=CONSTANT pi || (!isintriangle0 && !isintriangle2))
        && (isintriangle1 || isintriangle3) )
    {
        result[0] = 0;
        result[1] = 2;
        //if point lies in "left" triangle
        if(isintriangle1)
            result[2] = 3;
        //if point lies in "right" triangle
        else
            result[2] = 1;
    }
    //if the "left" and "right" angles sum to larger than 180
    else if(isintriangle0 || isintriangle2)
    {
        result[0] = 1;
        result[1] = 3;
        //if point lies in "top" triangle
        if(isintriangle0)
            result[2] = 0;
        //if point lies in "bottom" triangle
        else
            result[2] = 2;
    }
}
bool pixsrc_geometry::fallsinsrcinputcircle(double xx, double yy, double **circles, double size)
{
    for(PS_unsignedSIT x=0; x<size; x++)
        if( OPERA is_exactly_p_inf(circles[x][2]) )
            return 1;

    for(PS_unsignedSIT x=0; x<size; x++)
        if( ( (circles[x][0]-xx)*(circles[x][0]-xx) + (circles[x][1]-yy)*(circles[x][1]-yy) )
            <= circles[x][2] )
            return 1;

    return 0;
}
void pixsrc_geometry::getconvexhull( struct triangulateio *triinloc, struct triangulateio *trioutloc, PS_SIT *numvertsdest, double **convexhull_, char *switches )
{
    // any set of switches other than those in pixsrc_constant
    // are assumed to change the pointlist.
    // switches=NULL means switches defined in CONSTANT are used.
    // triinloc=NULL means no triangulation to be performed.
    // numvert=NULL means first element of convexhull will
    // contain number of vertices.

    if( !switches )
        switches = (char*)CONSTANT triswitchesnominangle;

    double *list;
    if( triinloc )
    {
        list = triinloc->pointlist;
        ps_tri_triangulate( switches, triinloc, trioutloc, (struct triangulateio*)NULL );
    }

    if( !triinloc || strcmp( switches, (char*)CONSTANT triswitchesnominangle) )
        list = trioutloc->pointlist;

    // get convex hull
    PS_SIT numverts = 0;
    for(PS_SIT j=0; j<trioutloc->numberofedges; j++)
        if(trioutloc->edgemarkerlist[j]==1)
        {
            ++numverts;
        }

    double *hullunsorted;
    MEMORY ps_malloc( &(hullunsorted), numverts*4 );
    PS_SIT index=0;
    for(PS_SIT j=0; j<trioutloc->numberofedges; j++)
    {
        if(trioutloc->edgemarkerlist[j]==1)
        {
            hullunsorted[index*2  ] = list[2*trioutloc->edgelist[2*j  ]  ];
            hullunsorted[index*2+1] = list[2*trioutloc->edgelist[2*j  ]+1];
            ++index;
            hullunsorted[index*2  ] = list[2*trioutloc->edgelist[2*j+1]  ];
            hullunsorted[index*2+1] = list[2*trioutloc->edgelist[2*j+1]+1];
            ++index;
        }
    }

    PS_SIT nelements = numverts*2;
    if( !numvertsdest )
        ++nelements;

    char *used;
    double *convexhull;
    MEMORY ps_malloc( &(*convexhull_), nelements );
    MEMORY ps_malloc( &(    used    ), numverts  );
    std::fill(used,used+numverts,0);


    if( !numvertsdest )
    {
        convexhull     = *convexhull_+1;
        convexhull[-1] = numverts;
    }
    else
    {
        convexhull    = *convexhull_;
        *numvertsdest = numverts;
    }

    // connects lines segments end to end, so that
    // convex hull is CW or CCW
    std::copy( hullunsorted, hullunsorted+4, convexhull );
    used[0]=1;
    for(PS_SIT k=2; k<numverts; k++)
    {
        for(PS_SIT j=0; j<2*numverts; j++)
        {
            if( !used[j/2] &&
                convexhull[(k-1)*2  ]==hullunsorted[j*2  ] &&
                convexhull[(k-1)*2+1]==hullunsorted[j*2+1] )
            {
                index = (j%2==0) ? j+1 : j-1;

                convexhull[k*2  ] = hullunsorted[index*2  ];
                convexhull[k*2+1] = hullunsorted[index*2+1];

                used[j/2]=1;
                break;
            }
        }
    }

    MEMORY ps_free( hullunsorted );
    MEMORY ps_free(     used     );
}

PS_SIT pixsrc_geometry::search_triangle( struct triangulateio* tri, PS_SIT *tri_seed, double x, double y )
{
    PS_SIT tri_i = (*tri_seed<0) ? 0 : *tri_seed;
    PS_SIT f_0, f_opp, iter=0;
    double det, t[2], ctr[2], v[2][3];
    while( 1 )
    {
        // if an infinite loop has occured while trying to find triangle that contains (x,y)
        if( iter>=tri->numberoftriangles )
            break;

        // setup triangle vertices
        for( PS_SIT f=0; f<3; ++f )
        {
            v[0][f] = tri->pointlist[tri->trianglelist[tri_i*3+f]*2  ];
            v[1][f] = tri->pointlist[tri->trianglelist[tri_i*3+f]*2+1];
        }

        // test if point lies in this triangle
        if( GEOM isintripadded( v, x, y ) )
        {
            *tri_seed = tri_i;
            return tri_i;
        }
        // compute next triangle to try
        else
        {
            // center of tri_i
            ctr[0] = ( v[0][0]+v[0][1]+v[0][2] ) / 3.0;
            ctr[1] = ( v[1][0]+v[1][1]+v[1][2] ) / 3.0;
            f_opp = 1;
            f_0   = 2;
            for( PS_SIT f=0; f<3; ++f )
            {
                // checking for edge of tri_i that intersects with line segment ctr--(x,y)
                // write the equations for the line segments parametrically:
                // line segment of ctr--(x,y) is given by: LS_ctr(t1)  = ctr + t1*((x,y)-ctr)
                // line segment of    edge    is given by: LS_edge(t2) = p1  + t2*( p2  -p1 )
                // where t1 and t2 are between 0 and 1, and p1 and p2 are endpoints
                // Then, after equating the two lines, and some algebra:

                det = (x-ctr[0])*(v[1][f]-v[1][f_0]) - (v[0][f]-v[0][f_0])*(y-ctr[1]);
                if( OPERA equalszero(det) )
                    continue;

                t[0] = ( (v[0][f_0]-ctr[0]) * (v[1][f]-v[1][f_0])
                         - (v[1][f_0]-ctr[1])*(v[0][f]-v[0][f_0]) ) / det;
                t[1] = ( (v[0][f_0]-ctr[0]) * (y-ctr[1])
                         - (v[1][f_0]-ctr[1])*(x-ctr[0]) ) / det;

                if( t[0]>=0 && t[0]<=1 && t[1]>=0 && t[1]<=1 )
                {
                    tri_i = tri->neighborlist[tri_i*3+f_opp];
                    if (tri_i!=-1)
                        break;
                    else return -1;
                }

                ++f_0;
                ++f_opp;
                if( f_0==3 )
                    f_0 = 0;
                if( f_opp==3 )
                    f_opp = 0;
            }
        }
        ++iter;
    }

    // if triangulation search failed for some reason, do brute force search
    for(PS_SIT tr = 0; tr < tri->numberoftriangles; ++tr)
    {
        for(PS_SIT f=0; f<3; ++f)
        {
            v[0][f] = tri->pointlist[tri->trianglelist[tr*3+f]*2  ];
            v[1][f] = tri->pointlist[tri->trianglelist[tr*3+f]*2+1];
        }
        if( GEOM isintripadded( v, x, y ) )
        {
            *tri_seed = tr;
            return tr;
        }
    }

    // if triangle still couldn't be found, point is outside trianuglation
    return -1;
}



// code from rosetta
typedef struct { double x, y; } vec_t;

inline double dot(vec_t *a, vec_t *b)
{
    return a->x * b->x + a->y * b->y;
}

inline double cross(vec_t *a, vec_t *b)
{
    return a->x * b->y - a->y * b->x;
}

inline vec_t *vsub(vec_t *a, vec_t *b, vec_t *res)
{
    res->x = a->x - b->x;
    res->y = a->y - b->y;
    return res;
}

/* tells if vec c lies on the left side of directed edge a->b
 * 1 if left, -1 if right, 0 if colinear
 */
PS_SIT left_of(vec_t *a, vec_t *b, vec_t *c)
{
    vec_t tmp1, tmp2;
    double x;
    vsub(b, a, &tmp1);
    vsub(c, b, &tmp2);
    x = cross(&tmp1, &tmp2);
    return x < 0 ? -1 : x > 0;
}

PS_SIT line_sect(vec_t *x0, vec_t *x1, vec_t *y0, vec_t *y1, vec_t *res)
{
    vec_t dx, dy, d;
    vsub(x1, x0, &dx);
    vsub(y1, y0, &dy);
    vsub(x0, y0, &d);
    /* x0 + a dx = y0 + b dy ->
       x0 X dx = y0 X dx + b dy X dx ->
       b = (x0 - y0) X dx / (dy X dx) */
    double dyx = cross(&dy, &dx);
    if (!dyx) return 0;
    dyx = cross(&d, &dx) / dyx;
    if (dyx <= 0 || dyx >= 1) return 0;

    res->x = y0->x + dyx * dy.x;
    res->y = y0->y + dyx * dy.y;
    return 1;
}

/* === polygon stuff === */
typedef struct { PS_SIT len, alloc; vec_t *v; } poly_t;

poly_t* poly_new()
{
    return (poly_t*)calloc(1, sizeof(poly_t));
}

void poly_free(poly_t *p)
{
    free(p->v);
    free(p);
}

void poly_append(poly_t *p, vec_t *v)
{
    if (p->len >= p->alloc) {
        p->alloc *= 2;
        if (!p->alloc) p->alloc = 4;
        p->v = (vec_t*)realloc(p->v, sizeof(vec_t) * p->alloc);
    }
    p->v[p->len++] = *v;
}

/* this works only if all of the following are true:
 *   1. poly has no colinear edges;
 *   2. poly has no duplicate vertices;
 *   3. poly has at least three vertices;
 *   4. poly is convex (implying 3).
 */
PS_SIT poly_winding(poly_t *p)
{
    return left_of(p->v, p->v + 1, p->v + 2);
}

void poly_edge_clip(poly_t *sub, vec_t *x0, vec_t *x1, PS_SIT left, poly_t *res)
{
    PS_SIT i, side0, side1;
    vec_t tmp;
    vec_t *v0 = sub->v + sub->len - 1, *v1;
    res->len = 0;

    side0 = left_of(x0, x1, v0);
    if (side0 != -left) poly_append(res, v0);

    for (i = 0; i < sub->len; i++) {
        v1 = sub->v + i;
        side1 = left_of(x0, x1, v1);
        if (side0 + side1 == 0 && side0)
            /* last point and current straddle the edge */
            if (line_sect(x0, x1, v0, v1, &tmp))
                poly_append(res, &tmp);
        if (i == sub->len - 1) break;
        if (side1 != -left) poly_append(res, v1);
        v0 = v1;
        side0 = side1;
    }
}

poly_t* poly_clip(poly_t *sub, poly_t *clip)
{
    PS_SIT i;
    poly_t *p1 = poly_new(), *p2 = poly_new(), *tmp;

    PS_SIT dir = poly_winding(clip);
    poly_edge_clip(sub, clip->v + clip->len - 1, clip->v, dir, p2);
    for (i = 0; i < clip->len - 1; i++) {
        tmp = p2; p2 = p1; p1 = tmp;
        if(p1->len == 0) {
            p2->len = 0;
            break;
        }
        poly_edge_clip(p1, clip->v + i, clip->v + i + 1, dir, p2);
    }

    poly_free(p1);
    return p2;
}

void pixsrc_geometry::clippoly (double* poly,     PS_SIT numpoly,
                                double* ps_clipper,  PS_SIT numclipper,
                                double** clipped, PS_SIT* numclipped)
{
    vec_t *c, *s;
    PS_SIT clen = numclipper;
    PS_SIT slen = numpoly;
    MEMORY ps_malloc (&c, numclipper);
    MEMORY ps_malloc (&s, numpoly);

    for (PS_SIT z=0; z<clen; ++z)
    {
        c[z].x = ps_clipper[z*2];
        c[z].y = ps_clipper[z*2+1];
    }
    for (PS_SIT z=0; z<slen; ++z)
    {
        s[z].x = poly[z*2];
        s[z].y = poly[z*2+1];
    }

    PS_SIT i;

    poly_t clipper = {clen, 0, c};
    poly_t subject = {slen, 0, s};

    poly_t *res = poly_clip(&subject, &clipper);

    *numclipped = res->len;
    MEMORY ps_malloc (&(*clipped), *numclipped*2);
    for (i = 0; i < res->len; i++)
    {
        (*clipped)[i*2]   = res->v[i].x;
        (*clipped)[i*2+1] = res->v[i].y;
    }

    poly_free (res);
    MEMORY ps_free (c);
    MEMORY ps_free (s);
}

void pixsrc_geometry::create_border (double *chull, PS_SIT numvert_chull,
                                     double spacing, double **border)
{
    // this function takes a convex hull and creates a list of points
    // that step along the hull, in steps of size "spacing"
    // resulting list of points is stored in "border"
    // the first element of border will hold the number of points

    // hacks
    PS_SIT i=0;
    PS_SIT numpts[1] = {0};
    double **polys = &chull;
    PS_SIT dim2 = numvert_chull*2;
    PS_SIT img_index = 1;
    double steps = spacing;

    // actual calculations starting

    for (PS_SIT k=0; k<dim2; k+=2)
    {
        PS_SIT j = (k==0) ? dim2-2 : k-2;

        // slope
        double m=0;
        if (OPERA equalszero(polys[i][k]-polys[i][j]))
            OPERA assign_p_infinity(&m);
        else
            m = (polys[i][k+1]-polys[i][j+1]) / (polys[i][k]-polys[i][j]);

        double coords[4] = {polys[i][k],polys[i][k+1],polys[i][j],polys[i][j+1]};
        // consecutive points that are exactly the same is a flag
        if (coords[0]==coords[2] && coords[1]==coords[3])
            continue;

        if (std::abs(m)<1)
        {
            if (polys[i][j]<polys[i][k])
            {
                coords[0] = polys[i][j];
                coords[1] = polys[i][j+1];
                coords[2] = polys[i][k];
                coords[3] = polys[i][k+1];
            }

            PS_SIT iter = 0;
            while (1)
            {
                // (x2-x1)^2 + (y2-y1)^2 = N^2*d^2
                // (x2-x1)^2 + m^2*(x2-x1)^2 = N^2*d^2
                // (1+m^2) * (x2-x1)^2 = N^2*d^2
                // x2 = [ N^2*d^2 / (1+m^2) ]^0.5 + x1
                double newx = coords[0] + iter*steps*std::sqrt(1.0 / (1.0+m*m));
                if (newx>=coords[2])
                    break;

                ++numpts[i];
                ++iter;
            }
        }
        else
        {
            if (polys[i][j+1]<polys[i][k+1])
            {
                coords[0] = polys[i][j];
                coords[1] = polys[i][j+1];
                coords[2] = polys[i][k];
                coords[3] = polys[i][k+1];
            }

            PS_SIT iter = 0;
            while (1)
            {
                // (x2-x1)^2 + (y2-y1)^2 = N^2*d^2
                // m^-2(y2-y1)^2 + (y2-y1)^2 = N^2*d^2
                // (1+m^-2) * (y2-y1)^2 = N^2*d^2
                // y2 = [ N^2*d^2 / (1+m^-2) ]^0.5 + y1
                double newy = coords[1] + iter*steps*std::sqrt(1.0 / (1.0+1.0/(m*m)));
                if (newy>=coords[3])
                    break;

                ++numpts[i];
                ++iter;
            }
        }
    }

    MEMORY ps_malloc (&(border[i]), numpts[i]*2+1);
    border[i][0] = numpts[i];

    // iterate over vertices
    // j is th "previous" vertex
    for (PS_SIT k=0; k<dim2; k+=2)
    {
        PS_SIT j = (k==0) ? dim2-2 : k-2;

        // slope
        double m=0;
        if (OPERA equalszero(polys[i][k]-polys[i][j]))
            OPERA assign_p_infinity(&m);
        else
            m = (polys[i][k+1]-polys[i][j+1]) / (polys[i][k]-polys[i][j]);
        double b = polys[i][k+1] - m*polys[i][k];

        double coords[4] = {polys[i][k],polys[i][k+1],polys[i][j],polys[i][j+1]};
        if (coords[0]==coords[2] && coords[1]==coords[3])
            continue;

        if (std::abs(m)<1)
        {
            if (polys[i][j]<polys[i][k])
            {
                coords[0] = polys[i][j];
                coords[1] = polys[i][j+1];
                coords[2] = polys[i][k];
                coords[3] = polys[i][k+1];
            }

            PS_SIT iter = 0;
            while (1)
            {
                double newx = coords[0] + iter*steps*std::sqrt(1.0 / (1.0+m*m));
                if (newx>=coords[2])
                    break;
                double newy = m*newx+b;

                border[i][img_index++] = newx;
                border[i][img_index++] = newy;

                ++iter;
            }
        }
        else
        {
            if (polys[i][j+1]<polys[i][k+1])
            {
                coords[0] = polys[i][j];
                coords[1] = polys[i][j+1];
                coords[2] = polys[i][k];
                coords[3] = polys[i][k+1];
            }

            PS_SIT iter = 0;
            while (1)
            {
                double newy = coords[1] + iter*steps*std::sqrt(1.0 / (1.0+1.0/(m*m)));
                if (newy>=coords[3])
                    break;
                double newx = (newy-b)/m;

                border[i][img_index++] = newx;
                border[i][img_index++] = newy;

                ++iter;
            }
        }
    }
}
