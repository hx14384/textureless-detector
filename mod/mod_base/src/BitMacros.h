/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

/**
*  @file BitMacros.h
*  @brief There are some macros for bit operations.
*/

#ifndef BITMACROS_H
#define BITMACROS_H
#define IS_SET_BIT(x,m)   ( ( (x) & (m) ) == (m) )
#define TOGGLE_BIT(x,m)   ( (x) ^=  (m) )
#define SET_BIT(x,m)      ( (x) |= (m) )
#define UNSET_BIT(x,m)    ( (x) &= (~m) )
#endif

