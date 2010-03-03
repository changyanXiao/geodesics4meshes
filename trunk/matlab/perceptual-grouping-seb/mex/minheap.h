//================================================================
//================================================================
// File: meanheap.h
// (C) 2008 by F. Benmansour and S. Bougleux
//================================================================
//================================================================
#ifndef __MINHEAP_H__
#define __MINHEAP_H__
//#include<limits>
//#include<iostream>
#include<cmath>
//================================================================
//================================================================
template<class Integer, class T>
class MinHeap
{
 private:
  Integer size;
  Integer *ptrToData;          // current position
  Integer *data;               // data vector
  T *_data;                    // data for classification rule
  Integer *tree;               // tree structure
  Integer Parent(Integer);
  Integer RightSon(Integer);
  Integer LeftSon(Integer);
  void _Update();
  //int TreePullFrom(int);
 public:
  MinHeap(unsigned long,T*);
  ~MinHeap();
  bool Empty();
  Integer Size();
  Integer Push(Integer);
  Integer Pop();
  void Update(Integer);
  void Clear();
  Integer operator[](Integer) const;
};

//================================================================
template<class Integer, class T> 
MinHeap<Integer,T>::MinHeap(unsigned long s, T *d)
  : size(s), ptrToData(0), data(new Integer[s]), _data(d), 
     tree(new Integer[s])
{
  ptrToData = data - 1;
}

//================================================================
template<class Integer, class T> 
MinHeap<Integer,T>::~MinHeap()
{
  if (tree) delete tree;
  if (data) delete data;
  _data = NULL;
}

//================================================================
template<class Integer, class T>
Integer MinHeap<Integer,T>::Parent(Integer pos)
{
  if (pos)
    {
      if (std::ceil((float)pos/2) == (float)pos/2)
	return (pos/2 - 1);
      return ((pos-1)/2);
    }
  return 0;
}

//================================================================
template<class Integer, class T>
Integer MinHeap<Integer,T>::Push(Integer new_pos)
{
  *(++ptrToData) = new_pos;
  Integer pos = (Integer)(ptrToData - data), pos1, pos2, par;
  tree[new_pos] = pos;
  Integer s_idx = data[pos], f_idx = data[Parent(pos)];
  while ((pos != 0) && (_data[s_idx] < _data[f_idx]))
    {
      pos1 = data[pos];
      par = Parent(pos);
      pos2 = data[par];
      tree[pos1] = par;
      data[pos] = pos2;
      tree[pos2] = pos;
      data[par] = pos1;
      pos = par;
      s_idx = data[pos];
      f_idx = data[Parent(pos)];
    }
  return (ptrToData - data + 1);
}

//================================================================
template<class Integer, class T>
bool MinHeap<Integer,T>::Empty()
{
  return ((ptrToData - data + 1) == 0);
}

//================================================================
template<class Integer, class T>
Integer MinHeap<Integer,T>::RightSon(Integer pos)
{
  if ((2 * pos + 2) < (Integer)(ptrToData - data + 1))
    return (2 * pos + 2);
  return 0;
}

//================================================================
template<class Integer, class T>
Integer MinHeap<Integer,T>::LeftSon(Integer pos)
{
  if ((2 * pos + 1) < (Integer)(ptrToData - data + 1))
    return (2 * pos + 1);
  return 0;
}

//================================================================
template<class Integer, class T>
void MinHeap<Integer,T>::_Update()
{
  Integer position = 0, ls_idx, rs_idx, idx, idx_r, idx_l, par;
  bool stop = false;
  while ((position >= 0) && (!stop))
    {		
      if(((idx_r = RightSon(position)) > 0) && 
	 ((idx_l = LeftSon(position)) > 0))
	{
	  ls_idx = data[idx_l];
	  rs_idx = data[idx_r];
	  if(_data[ls_idx] <= _data[rs_idx])
	    {
	      idx = data[position] = ls_idx;
	      tree[idx] = position;
	      position = idx_l;
	    }
	  else
	    {
	      idx = data[position] = rs_idx;
	      tree[idx] = position;
	      position = idx_r;
	    }
	}
      else
	{
	  if ((idx_l = LeftSon(position)) > 0)
	    {
	      idx = data[position] = data[idx_l];
	      tree[idx] = position;
	      position = idx_l;
	    }
	  else stop = true;
	}
    }
    if (position != (ptrToData - data))
      {
	tree[*ptrToData] = position;
        ls_idx = data[position] = *ptrToData;
        rs_idx = data[Parent(position)];
        while ((position != 0) && (_data[ls_idx] < _data[rs_idx]))
	  {
	    idx = data[position];				
	    par = tree[idx] = Parent(position);
	    ls_idx = data[position] = data[par];	
	    tree[ls_idx] = position;
	    data[par] = idx;
	    position = par;
	    rs_idx = data[Parent(position)];
	  }
      }
}

//================================================================
template<class Integer, class T>
Integer MinHeap<Integer,T>::Pop()
{
  if (ptrToData - data + 1)
    {
      Integer first = *data;
      tree[first] = 0;
      _Update();
      ptrToData--;
      return first;
    }
  return 0;
}

//================================================================
template<class Integer, class T>
void MinHeap<Integer,T>::Update(Integer position)
{
  Integer s_idx = data[position], 
    f_idx = data[Parent(position)], idx, par;
  while ((position != 0) && (_data[s_idx] < _data[f_idx]))
    {
      idx = data[position];
      par = tree[idx] = Parent(position);
      s_idx = data[position] = data[par];
      tree[s_idx] = position;
      data[par] = idx;
      position = par;
      f_idx = data[Parent(position)];
    }
}

//================================================================
//template<class Integer, class T>
//int Tree::PullFrom(int point)
//================================================================
//{
//T Uv = _data[point];
//_data[point] = 0;
//Update(tree[point]);
//_data[Pop()] = Uv;
//}

//================================================================
template<class Integer, class T>
Integer MinHeap<Integer,T>::Size()
{
  return (ptrToData - data + 1);
}

//================================================================
template<class Integer, class T>
void MinHeap<Integer,T>::Clear()
{
  while (!Empty())
    {
      Pop();
    }
  //for (Integer i = 0; i < size; i++) data[i] = 0;
  ptrToData = data - 1;
}

//================================================================
template<class Integer, class T>
Integer MinHeap<Integer,T>::operator[](Integer pos) const
{
  return tree[pos];
}

#endif
