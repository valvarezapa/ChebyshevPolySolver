# Last update: 06/10/23 ---------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#----------------------------------- COMBINATORICS MODULE -----------------------------------------#
#--------------------------------------------------------------------------------------------------#
#-------- This module contains algorithms of combinatorial nature able to undertake ---------------#
#-------- calculations with cycles, compute inverses, discern whether a cycle  --------------------#
#-------- is in an specific conjugacy class,... along with some other methods that ----------------#
#-------- deal with lists, of general utility. ----------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

module com 

using Colors
using PyPlot


# Computes a vector containing the positions in which 'element' appears in the given list 'list'.
function positions(list,element)
    # Pre:
    # Post: Returns the vector of positions.
    positions_vector=[]
    l=length(list)
    @inbounds for i in 1:l
      if list[i]==element
        push!(positions_vector,i)
      end
    end
    return positions_vector
  end
  
  
  #Example
  #positions([1,2,4,6,1,25,3,8,4,2,4],4)


# Computes a vector containing the positions in which 'element' appears in the given list 'list'.
function positions(list::Array{Complex{Float64},1},element::Complex{Float64};prec::Float64=8.0)
  # Pre:
  # Post: Returns the vector of positions.
  positions_vector=[]
  l=length(list)
  @inbounds for i in 1:l
    if abs(list[i]-element)<1.0/10.0^prec
      push!(positions_vector,i)
    end
  end
  return positions_vector
end


#Example
#positions([1,2,4,6,1,25,3,8,4,2,4],4)


# 
function is_in(el::Complex{Float64},list::Array{Complex{Float64},1};pre::Float64=8.0)::Bool
    # Pre:
    # Post:
    l=length(list)
    for i in 1:l
      if abs(el-list[i])<1.0/(10.0^pre)
        return true
      end
    end
    return false
  end
  
  
  #Example
  #is_in()


  #
  function equal_list(list1::Array{Complex{Float64},1},list2::Array{Complex{Float64},1};precision::Float64=8.0)::Bool
    # Pre:
    # Post:
    l1=length(list1)
    l2=length(list2)
    if l1!=l2
      return false
    else
      for i in 1:l1
        if is_in(list1[i],list2;pre=precision)==false
          return false
        end
      end
      return true
    end
  end
  
  
  #Example
  #equal_list()


  # 
function delete_repited(list::Array{ComplexF64,1};precision::Float64=8.0)::Array{Complex{Float64},1}
    # Pre:
    # Post:
    l=length(list)
    new_list=[first(list)]
    for i in 2:l
      if is_in(list[i],new_list;pre=precision)==false
        push!(new_list,list[i])
      end
    end
    return new_list
  end
  
  
  #Example
  #delete_repited()


  # AVOIDS UNDERFLOWS
function clean_list(list::Array{Complex{Float64},1};precision::Float64=8.0)::Array{Complex{Float64},1}
    # Pre:
    # Post:
    l=length(list)
    new_list=[]
    for i in 1:l
      x=real(list[i])
      y=imag(list[i])
      if abs(x)<1.0/(10.0^precision)
        x=0
      end
      if abs(y)<1.0/(10.0^precision)
        y=0
      end
      push!(new_list,x+y*im)
    end
    return new_list
  end
  
  
  #Example
  #clean_list()

  #--------------------------------- CYCLES ----------------------------------------------#

  # 
function sameLengthCycles(cycle1,cycle2)
    # Pre: 
    # Post: 
    l1=length(cycle1)
    l2=length(cycle2)
    if l1!=l2
        if l1<l2
         while l1<l2
             l1=l1+1
             push!(list1,l1)
            end
        else
            while l2<l1
               l2=l2+1
               push!(list2,l2)
            end
        end
    end
  end
  
  
  # Computes the product of two cycles in the symmetric group S_n of order n.
  function cycleProduct(cycle1::Array{Int64,1},cycle2::Array{Int64,1})::Array{Int64,1}
    # Pre:
    # Post:
    sameLengthCycles(cycle1,cycle2)
    l=length(cycle1)
    result=[]
    for i in 1:l
      push!(result,cycle1[cycle2[i]])
    end
    return result
  end
  
  
  #Example
  #cycleProduct()
  
  
  # 
  function getPosition(el::Int64,cycle::Array{Int64,1})::Int64
    # Pre:
    # Post:
    l=length(cycle)
    for i in 1:l
      if cycle[i]==el
        return i
      end
    end
  end


  #Example
  #cycleProduct()
  
  
  # Computes the inverse element of a cycle in the symmetric group S_n of order n.
  function inverseCycle(cycle::Array{Int64,1})::Array{Int64,1}
    # Pre:
    # Post:
    l=length(cycle)
    inverse=[]
    for i in 1:l
      push!(inverse,getPosition(i,cycle))
    end
    return inverse
  end


  #Example
  #cycleProduct()


  # 
  function equalCycles(cycle1::Array{Int64,1},cycle2::Array{Int64,1})::Bool
    # Pre: The cycles have the same length.
    # Post:
    l=length(cycle1)
    for i in 1:l
      if cycle1[i]!=cycle2[i]
        return false
      end
    end
    return true
  end


  #Example
  #cycleProduct()


  # Computes whether if the two given cycles are conjugate in the symmetric group S_n of order n.
  function are_conj(cycle1::Array{Int64,1},cycle2::Array{Int64,1})::Array{Any,1} # ONLY FOR S3
  # Pre:
  # Post:
  sameLengthCycles(cycle1,cycle2)
  if length(cycle1)!=3
    return [[],false]
  else
    els=[]
    elements=[[1,2,3],[1,3,2],[2,1,3],[3,2,1],[2,3,1],[3,1,2]]
    for el in elements
      if equalCycles(cycleProduct(cycleProduct(el,cycle1),inverseCycle(el)),cycle2)
        push!(els,el)
      end
    end
    if length(els)==0
      return [els,false]
    else
      return [els,true]
    end
  end
  end


  #Example
  #are_conj([2,1,3],[1,3,2])


  # Computes whether if the two given elements in S_n x...(v)...x S_n are conjugate.
  function are_conjCanon(canon1::Array{Array{Int64,1},1},canon2::Array{Array{Int64,1},1})::Bool # ONLY FOR S3
    # Pre:
    # Post:
    l1=length(canon1)
    if l1!=length(canon2)
      return false
    end
    for i in 1:l1
      if length(canon1[i])!=length(canon2[i])
        return false
      end
    end
    conjFirstGens=are_conj(canon1[1],canon2[1])
    if conjFirstGens[2]==true
      elements=conjFirstGens[1]
    else
      return false
    end
    for el in elements
      areConj=true
      i=2
      while i<=l1 && areConj
        if equalCycles(com.cycleProduct(com.cycleProduct(el,canon1[i]),com.inverseCycle(el)),canon2[i])
          i=i+1
        else
          areConj=false
        end
      end
      if areConj==true
        return true
      end
    end
    return false
    end
  
  
    #Example
    #are_conj([2,1,3],[1,3,2])


  # Computes all the transitive conjugacy classes in  S_n x S_n, where S_n denotes the symmetric group of order n.
  function transConClasses(n::Int64)::Array{Any,1} # ONLY FOR S3
    # Pre:
    # Post:
    @assert n==3
    els=[[1,2,3],[1,3,2],[2,1,3],[3,2,1],[2,3,1],[3,1,2]]
    elements=[]
    for el1 in els
      for el2 in els
        push!(elements,[el1,el2])
      end
    end
    #print("Elements in S3xS3 (");print(length(elements));println("):")
    #println(elements)

    transels=[]
    for el in elements
      istrans=true
      i=1
      while i<=n && istrans==true
        if el[1][i]==i && el[2][i]==i
          istrans=false
        end
        i=i+1
      end
      if istrans==true
        push!(transels,el)
      end
    end
    #print("Transitive elements in S3xS3 (");print(length(transels));println("):")
    #println(transels)

    classes=[pop!(transels)]
    for el in transels
      classified=false
      num_classes=length(classes)
      i=1
      while i<=num_classes && classified==false
        areConjEls=are_conj(el[1],classes[i][1])
        areConj=areConjEls[2]
        num_pos=length(areConjEls[1])
        j=1
        while j<=num_pos && areConj && classified==false
          if equalCycles(com.cycleProduct(com.cycleProduct(areConjEls[1][j],el[2]),com.inverseCycle(areConjEls[1][j])),classes[i][2])
            classified=true
          else
            j=j+1
          end
        end
        i=i+1
      end
      if classified==false
        push!(classes,el)
      end
    end
    #print("Conjugacy classes of transitive elements in S3xS3 (");print(length(classes));println("):")
    #println(classes)

    return classes
  end
  
  
    #Example
    #are_conj([2,1,3],[1,3,2])


end