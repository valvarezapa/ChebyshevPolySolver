# Last update: 09/10/23 ---------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#----------------------------------- GEOMETRY MODULE ----------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#-------- This module contains algorithms of geometric nature able to construct -------------------#
#-------- paths on the Riemann sphere and to operate with them, aswell as to ----------------------#
#-------- tackle some other purely geometrical aspects of our model. ------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

module geo

include("combinatorics-lift.jl")

using Colors
using PyPlot


# Computes the norm of the maximum of a complex number.
function norm_max(z::Complex{Float64})::Float64
  # Pre:
  # Post: Returns the maximum between the real part and the imaginary part of the given complex number.
  return max(abs(real(z)),abs(imag(z)))
end
  
  
#Example
#norm_max(1.0-3.0*im)


# Computes the distance (considering the norm of the maximum) between two complex numbers.
function dist_max(z1::Complex{Float64},z2::Complex{Float64})::Float64
  # Pre:
  # Post: Returns the distance between two given complex numbers.
  return norm_max(z2-z1)
end


#Example
#dist_max(1.0-3.0*im,complex(25.0))


# Computes the distance (considering the norm of the maximum) between two complex numbers.
function dist_max(points::Tuple{Complex{Float64},Complex{Float64}})::Float64
  # Pre:
  # Post: Returns the distance between two given complex numbers.
  return norm_max(points[2]-points[1])
end


#Example
#dist_max((1.0-3.0*im,complex(25.0)))


# Computes the distance (considering the norm of the maximum) between a complex number and a list of points.
function dist_max(z::Complex{Float64},points::Array{Complex{Float64},1})::Float64
  # Pre:
  # Post: Returns the distance between the given point and the list of complex numbers.
  l=length(points)
  if l>0
    dis=dist_max(z,points[1])
    @inbounds for i in 2:l
      if dist_max(z,points[i])<dis
        dis=dist_max(z,points[i])
      end
    end
  else
    dis=0
  end
  return dis
end


#Example
#dist_max((1.0-3.0*im,complex(25.0)))


# Computes the minimum distance (considering the norm of the maximum) in a list of complex numbers.
function min_dist_max(points::Array{Complex{Float64},1})::Float64
  # Pre:
  # Post: Returns the minimum distance between points of the given list.
  l=length(points)
  if l>1
    min=dist_max(points[1],points[2])
    @inbounds for i in 1:l
      @inbounds for j in 1:l
        if i!=j && dist_max(points[i],points[j])<min
          min=dist_max(points[i],points[j])
        end
      end
    end
  else
    min=0
  end
  return min
end


#Example
#min_dist_max((1.0-3.0*im,complex(25.0),3.05+18*im))


# Computes the direction of a discrete path joining two complex numbers (z1 and z2), in two steps.
# The first component of the resulting vector is the direction which z1 has to follow in order to get to the intermediate point m1.
# The second component of the resulting vector is the direction which m1 has to follow in order to get to the next point m2 of the discrete path.
# In case the two points are vertically or horizontally aligned, the resulting vector only has one element.
# This choice of steps will determine the shape of the discrete paths we will be able to construct.
function direction(z1::Complex{Float64},z2::Complex{Float64})::Array{Complex{Float64},1}
  # Pre:
  # Post: Returns the direction vector.
  r1=real(z1)
  i1=imag(z1)
  r2=real(z2)
  i2=imag(z2)
  if abs(r2-r1)<abs(i2-i1) # If the horizontal difference between the two points is lesser than the vertical difference, the first step taken is vertical.
    if r2!=r1 
      return [((i2-i1)/abs(i2-i1))*im,complex((r2-r1)/abs(r2-r1))]
    else # In case the two points have the same real part (vertically aligned), the direction will be entirely vertical.
      return [((i2-i1)/abs(i2-i1))*im]
    end
  else # If the vertical difference between the two points is lesser or equal than the horizontal difference, the first step taken is horizontal.
    if i2!=i1
      return [complex((r2-r1)/abs(r2-r1)),((i2-i1)/abs(i2-i1))*im]
    else # In case the two points have the same imaginary part (horizontally aligned), the direction will be entirely horizontal.
      return [complex((r2-r1)/abs(r2-r1))]
    end
  end
end


#Example
#direction(1.0-3.0*im,complex(25.0))


# Constructs a discrete corner joining two complex numbers.
function corner(z1::Complex{Float64},z2::Complex{Float64})::Array{Complex{Float64},1}
  # Pre:
  # Post: Returns the list of points that make up the corner.
  r1=real(z1)
  i1=imag(z1)
  r2=real(z2)
  i2=imag(z2)
  if abs(r2-r1)<abs(i2-i1) # Depending on the position of both points we construct a suitable discrete corner joining them. 
    if r2!=r1 
      return [z1,r1+i2*im,z2]
    else 
      return [z1,z2]
    end
  else
    if i2!=i1
      return [z1,r2+i1*im,z2]
    else 
      return [z1,z2]
    end
  end
end


#Example
#corner(1.0-3.0*im,complex(25.0))


# Constructs a discrete path joining the two given complex numbers z1 and z2, avoiding to intersect an square-shaped neighborhood of each point in the given list of points.
function polygonal_path(points::Array{Complex{Float64},1},z1::Complex{Float64},z2::Complex{Float64})::Array{Complex{Float64},1}
  # Pre:
  # Post: Returns the list of points that make up the polygonal path.
  if z1==z2 # If the two points are the same, we construct a discrete loop based on z1 that intersects z1+1.
    pol_path=vcat(polygonal_path(points,z1,z1+1),polygonal_path(points,z1+1,z1))
    return pol_path
  else
    l=length(points)
    if l==0 # If there are no points to avoid, we simply construct a corner between z1 and z2. 
      return corner(z1,z2)
    else
      if l==1
        mult=minimum([dist_max(z1,points[1]),dist_max(z2,points[1])])/3.0 # The value of mult will determine the size of the steps in the constructed discrete path. Note that a sufficiently low value of mult
      else                                                                # will result in the discrete path having too many points, which might lead to some high execution time or even some numerical problems.
        mult=minimum([dist_max(z1,points),dist_max(z2,points),min_dist_max(points)])/3.0
      end
      if dist_max(z1,z2)<=mult # If z1 and z2 are closer to each other rather than to the points to avoid, we simply return a corner between them.
        return corner(z1,z2)
      else
        pol_path=[z1]
        dir=direction(last(pol_path),z2)
        while (dist_max(last(pol_path),z2)>mult) # While the last point in our path is far enough from the ending point z2.
          if (dist_max(last(pol_path),points)>=mult) # If the last point in our path is far enough from every point in the list, we can take a step in the direction of z2.
            dir=direction(last(pol_path),z2)
            @inbounds for i in 1:length(dir)
              push!(pol_path,last(pol_path)+mult*dir[i])
            end
          else # If the last point in our path is too close to some point in the list (in other words, it is in an square-shaped neighborhood of some point in the list).
            pop!(pol_path) # We take a step back and two steps forward in the same direction in order to avoid the point.
            if length(dir)==1
              dir=[dir[1]*im]
            else
              dir=[dir[1]]
            end
            push!(pol_path,last(pol_path)+1.1*mult*dir[1])
          end
        end
        cor=corner(last(pol_path),z2) # Finally, when the last point in our path is sufficiently close to z2, we make a corner between those points.
        len_cor=length(cor)
        @inbounds for i in 1:len_cor
          push!(pol_path,cor[i])
        end
        return pol_path
      end
    end
  end
end


#Example
#polygonal_path([complex(10.0),15.0+3*im],1.0-3.0*im,complex(25.0))


# Constructs and plots a discrete path joining the two given complex numbers, avoiding to intersect a neighborhood of each point in the given list.
function plot_polygonal_path(points::Array{Complex{Float64},1},z1::Complex{Float64},z2::Complex{Float64};plotcolor::AbstractString="blue")
  # Pre:
  # Post: Plots the polygonal path.
  len_points=length(points)
  pol_path=polygonal_path(points,z1,z2)
  len_path=length(pol_path) # An extremely large number of points in the polygonal path might lead to computational problems.
  if len_points==0 # If there are no points to avoid, we simply plot each point in the discrete path in the chosen color.
    @inbounds for i in 1:len_path-1
      x1=real(pol_path[i])
      y1=imag(pol_path[i])
      x2=real(pol_path[i+1])
      y2=imag(pol_path[i+1])
      x_values=[x1,x2]
      y_values=[y1,y2]
      PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
    end
  else # In other case, we will also construct and plot an square-shaped neighborhood of each point to avoid in gray color.
    if len_points==1
      mult=minimum([1.0,dist_max(z1,points[1]),dist_max(z2,points[1])])
    else
      mult=minimum([dist_max(z1,points),dist_max(z2,points),min_dist_max(points)/3.0])
    end
    @inbounds squares=[((real(points[k])-mult,imag(points[k])-mult),(real(points[k])+mult,imag(points[k])-mult),(real(points[k])+mult,imag(points[k])+mult),(real(points[k])-mult,imag(points[k])+mult)) for k=1:len_points]
    @inbounds for i in 1:len_points
      @inbounds for j in 1:4
        if j==4
          x1=squares[i][j][1]
          y1=squares[i][j][2]
          x2=squares[i][1][1]
          y2=squares[i][1][2]
        else
          x1=squares[i][j][1]
          y1=squares[i][j][2]
          x2=squares[i][j+1][1]
          y2=squares[i][j+1][2]
        end
        x_values=[x1,x2]
        y_values=[y1,y2]
        PyPlot.plot(x_values,y_values,color="gray",linestyle="-")
      end
    end
    @inbounds for i in 1:len_path-1
      x1=real(pol_path[i])
      y1=imag(pol_path[i])
      x2=real(pol_path[i+1])
      y2=imag(pol_path[i+1])
      x_values=[x1,x2]
      y_values=[y1,y2]
      PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
    end
  end
end


#Example
#clf() # Can be useful to clear the current plot before making a new one.

#plot_polygonal_path([complex(10.0)],1.0-3.0*im,complex(25.0))
#gcf()


# Plots a given discrete path, previously constructed.
function plot_polygonal_path(pol_path::Array{Complex{Float64},1};plotcolor::AbstractString="blue")
  # Pre:
  # Post: Plots the polygonal path.
  len_path=length(pol_path)
  @inbounds for i in 1:len_path-1
    x1=real(pol_path[i])
    y1=imag(pol_path[i])
    x2=real(pol_path[i+1])
    y2=imag(pol_path[i+1])
    x_values=[x1,x2]
    y_values=[y1,y2]
    PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
  end
end


# Constructs a discrete straight segment joining two complex numbers.
# If the value of the variable 'string' is 'c', the whole segment will be returned as an array of complex points.
# If 'string' is 'r', the last point of the segment (z2) will be removed, whereas if 'string' is 'l', the first
# point (z1) will be removed.
# If the value of the variable 'reverse' is -1, the program returns the corresponding reverse segment between z2 and z1.
function discrete_segment(z1::Complex{Float64},z2::Complex{Float64},num_points::Int64=20;string::Char='c',reverse::Int64=1)::Array{Complex{Float64},1}
  # Pre: 'reverse'={1,-1}, 'string'={'c','r','l'}
  # Post: Returns the list of points that make up the segment.
  if reverse==1
    @inbounds seg=[z1+(i/num_points)*(z2-z1) for i in 0:num_points]
    if string=='c'
      return seg
    elseif string=='r'
      pop!(seg)
      return seg
    elseif string=='l'
      popfirst!(seg)
      return seg
    end
  elseif reverse==-1
    if string=='c'
      return discrete_segment(z2,z1,num_points;string='c',reverse=1)
    elseif string=='r' 
      return discrete_segment(z2,z1,num_points;string='l',reverse=1)
    elseif string=='l'
      return discrete_segment(z2,z1,num_points;string='r',reverse=1)
    end
  end
end


#Example
#discrete_segment(1.0-3.0*im,complex(25.0),50;string='r',reverse=1)


# Constructs a discrete square loop of given 'radius' around a complex number.
# If the value of the variable 'string' is 'c', the whole loop will be returned as an array of complex points.
# If 'string' is 'r', the last point of the loop (center) will be removed, whereas if 'string' is 'l', the first
# point (center) will be removed.
# If the value of the variable 'reverse' is -1, the program returns the corresponding reversed loop.
function discrete_loop(center::Complex{Float64},radius::Float64,num_points::Int64=20;string::Char='c',reverse::Int64=1)::Array{Complex{Float64},1}
  # Pre: 'reverse'={1,-1}, 'string'={'c','r','l'}
  # Post: Returns the list of points that make up the loop.
  first_half_right=discrete_segment(center+complex(radius), center+complex(radius)+radius*im, num_points;string='r',reverse=1)
  up=discrete_segment(center+complex(radius)+radius*im, center-complex(radius)+radius*im, 2*num_points;string='r',reverse=1)
  left=discrete_segment(center-complex(radius)+radius*im, center-complex(radius)-radius*im, 2*num_points;string='r',reverse=1)
  down=discrete_segment(center-complex(radius)-radius*im, center+complex(radius)-radius*im, 2*num_points;string='r',reverse=1)
  second_half_right=discrete_segment(center+complex(radius)-radius*im, center+complex(radius), num_points;string='c',reverse=1)
  loop=[first_half_right; up; left; down; second_half_right]
  if reverse==1
    if string=='c'
      return loop
    elseif string=='r'
      pop!(loop)
      return loop
    elseif string=='l'
      popfirst!(loop)
      return loop
    elseif string=='o'
      popfirst!(loop)
      pop!(loop)
      return loop
    end
  elseif reverse==-1
    if string=='c'
      return reverse!(loop)
    elseif string=='r' 
      pop!(loop)
      return reverse!(loop)
    elseif string=='l'
      popfirst!(loop)
      return reverse!(loop)
    elseif string=='o'
      popfirst!(loop)
      pop!(loop)
      return reverse!(loop)
    end
  end
end


#Example
#discrete_loop(1.0-3.0*im,2.5,50;string='r',reverse=1)


# Computes the 'base point' of the given list of points.
# The base point is considered to be the real point with value equal to the maximum real part of the points in the list plus the minimum distance 
# (using the norm of the maximum) between the points.
function basepoint(points::Array{Complex{Float64},1})::Complex{Float64}
  # Pre:
  # Post: Returns the base point.
  l=length(points)
  if l==1
    return complex(real(points[1])+1.0/3.0)
  elseif l==0
    return complex(0.0)
  else
    @inbounds points_real=[real(points[i]) for i in 1:l]
    xmax=maximum(points_real)
    xmax1=xmax+min_dist_max(points)
    return complex(xmax1)
  end
end


#Example
#basepoint([1.0-3.0*im,complex(2.5),50.0+6.2*im])


# Computes a discrete loop that starts and ends at the base point of the given list of points, and that loops around the 'loop_center',
# avoiding to intersect an square-shaped neighborhood of radius 'radius' of the points in the given list.
function discrete_loop(points::Array{Complex{Float64},1},loop_center::Complex{Float64},radius::Float64,num_points::Int64=20)::Array{Complex{Float64},1}
  # Pre:
  # Post: Returns the discrete loop starting and ending at 'loop_center'.
  bp=basepoint(points)
  segment=polygonal_path(points, bp, loop_center+complex(radius))
  loop_around=discrete_loop(loop_center,radius,num_points;string='o',reverse=1)
  return [segment; loop_around; reverse(segment)]
end


#Example
#discrete_loop([1.0-3.0*im,complex(2.5),50.0+6.2*im],complex(0.0),0.5)


# Computes a list of discrete loops that start and end at the base point of the given list of points, and that loops around the each point of the list, 
# avoiding to intersect an square-shaped neighborhood of the  rest of the points in the given list.
function discrete_loops(points::Array{Complex{Float64},1},num_points::Int64=20)::Array{Array{Complex{Float64},1},1}
  # Pre:
  # Post: Returns the list of discrete loops.
  loops=[]
  l=length(points)
  if l==1
    radius=1.0/3.0
  elseif l==0
    return loops
  else
    radius=min_dist_max(points)/3.0
  end
  @inbounds for i in 1:l
    new_loop=discrete_loop(points,points[i],radius,num_points)
    push!(loops,new_loop)
  end
  return loops
end


#Example
#discrete_loops([1.0-3.0*im,complex(2.5),50.0+6.2*im])


# Plots a given loop, avoiding to intersect an square-shaped neighborhood of each point in the given list.
function plot_discrete_loop(points::Array{Complex{Float64},1},loop::Array{Complex{Float64},1};plotcolor::AbstractString="blue",enable_plot::Bool=true)::Array{Any,1}
  # Pre: The loop does not intersect a fixed neighborhood of any of the points in the given list.
  # Post: Plots the loop and the obstacles, and returns their elements joined in a list.
  l=length(points)
  len_loop=length(loop)
  if l==0
    if enable_plot==true
      @inbounds for i in 1:len_loop-1
        x1=real(loop[i])
        y1=imag(loop[i])
        x2=real(loop[i+1])
        y2=imag(loop[i+1])
        x_values=[x1,x2]
        y_values=[y1,y2]
        PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
      end
    end
    return loop
  else
    mult=min_dist_max(points)/3.0
    @inbounds squares=[((real(points[k])-mult,imag(points[k])-mult),(real(points[k])+mult,imag(points[k])-mult),(real(points[k])+mult,imag(points[k])+mult),(real(points[k])-mult,imag(points[k])+mult)) for k=1:l]
    if enable_plot==true
      @inbounds for i in 1:l
        @inbounds for j in 1:4
          if j==4
            x1=squares[i][j][1]
            y1=squares[i][j][2]
            x2=squares[i][1][1]
            y2=squares[i][1][2]
          else
            x1=squares[i][j][1]
            y1=squares[i][j][2]
            x2=squares[i][j+1][1]
            y2=squares[i][j+1][2]
          end
          x_values=[x1,x2]
          y_values=[y1,y2]
          PyPlot.plot(x_values,y_values,color="gray",linestyle="-")
        end
      end
      @inbounds for i in 1:len_loop-1
        x1=real(loop[i])
        y1=imag(loop[i])
        x2=real(loop[i+1])
        y2=imag(loop[i+1])
        x_values=[x1,x2]
        y_values=[y1,y2]
        PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
      end
    end
    return [loop; squares]
  end
end


#Example
#plot_discrete_loop([1.0-3.0*im,complex(2.5),50.0+6.2*im],loop)


# Given a set of loops based on the same point and a path that ends in that base point, the algorithm computes a set of loops
# pair-wise homotopic to the given set of loops based on the starting point of the given path.
function base_change(path::Array{Complex{Float64},1},loops::Array{Array{Complex{Float64},1},1})::Array{Array{Complex{Float64},1},1}
  # Pre: Every loop in 'loops' is based on the same point, which is the end of the given path 'path'. 
  # Post: Returns the new set of loops, based on the new base point (the starting point of 'path').
  len_loops=length(loops)
  len_path=length(path)
  @inbounds path_aux=[path[k] for k in 1:len_path]
  pop!(path_aux)
  new_loops=[]
  @inbounds for i in 1:len_loops
    push!(new_loops,[path_aux; loops[i]; reverse(path_aux)])
  end
  return new_loops
end


#Example
#base_change()


# 
function non_quasicritical_value(value::Complex{Float64},critical_values::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Any,1}
  # Pre:
  # Post:
  l=length(critical_values)
  @inbounds distances=[dist_max(value,critical_values[k]) for k in 1:l]
  if l==0
    return [true] # 'value' is not among the given critical values.
  else
    min=minimum(distances)
    if min<1.0/10.0^precision
      positions_list=com.positions(distances,min)
      critical_values_list=[false,positions_list[1],critical_values[positions_list[1]]] # 'value' is the critical value in the position positions_list[1].
      return critical_values_list
    else
      return [true] # 'value' is not among the given critical values.
    end
  end
end


#Example
#non_quasicritical_value()


# Builds a collection of points of a rectangle determined by the given intervals and precision,
# respecting the aspect ratio between the given intervals. The given precision is used to
# determine the number of points in the grid.
# Note that the collection contains P^1+(C) points [z:t] that, considering the P^1+(C) to C*
# (that is the Alexandroff compactification of C) bijection, correspond to z in C*.
function rectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
  precision=3)
  # Pre:
  # Post: Returns the collection of points.
  resolution=9.443329*(10.0^(precision)) # We take this 9-million-pixels-resolution as a reference.
  a=xinterval[1]
  b=xinterval[2]
  c=yinterval[1]
  d=yinterval[2]
  p1=floor(sqrt(resolution*((b-a)/(d-c)))) # Note that resolution=p1*p2
  p2=floor(p1*(d-c)/(b-a))
  collection=[complex(r,i) for i=d:-(d-c)/p2:c, r=a:(b-a)/p1:b]
  return collection
end


#Example
#rectangle()


# Establishes a canonical fixed bijection between a set of n complex values and the canonical group of n elements {1,2,...,n}.
function fixedBijection(points::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Complex{Float64},1}
  # Pre:
  # Post:
  ang_coord=[]
  new_points=[]
  for p in points
    x=real(p)
    y=imag(p)
    if x>0 && y>=0
      push!(ang_coord,atan(y/x))
    elseif x==0 && y>0
      push!(ang_coord,pi/2)
    elseif x<0
      push!(ang_coord,atan(y/x)+pi)
    elseif x==0 && y<0
      push!(ang_coord,3*pi/2)
    else
      push!(ang_coord,atan(y/x)+2*pi)
    end
  end
  while length(points)>0
    min=minimum(ang_coord)
    pos=com.positions(ang_coord,min)
    if length(pos)>1
      rad_coord=[]
      for i in pos
        push!(rad_coord,sqrt(real(points[i])^2+imag(points[i])^2))
      end
      while length(rad_coord)>1
        minRad=minimum(rad_coord)
        posRad=com.positions(rad_coord,minRad)
        push!(new_points,points[pos[posRad[1]]])
        deleteat!(rad_coord,posRad[1])
        deleteat!(points,pos[posRad[1]])
      end
    elseif length(pos)==1
      push!(new_points,points[pos[1]])
      deleteat!(ang_coord,pos[1])
      deleteat!(points,pos[1])
    else
      println("Error in method fixedBijection")
    end
  end
  return new_points
end


#Example
#non_quasicritical_value()

end