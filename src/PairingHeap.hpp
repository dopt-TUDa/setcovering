// -*- C++ -*-
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    Implementation of pairing heaps                                        */
/*                                                                           */
/*    Copyright (C) 2003  Marc Pfetsch                                       */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PAIRING_HEAP_H
#define PAIRING_HEAP_H


/**@file   PairingHeap.hpp
 * @brief  Provides a pairing heap
 * @author Marc Pfetsch
 *
 * Pairing heaps are a fast alternative to Fibonacci heaps. A pairing heap is easier
 * to implement and much faster in practice than a Fibonacci heap. However, the asymptotic
 * complexity of some most operations is not known exactly.
 *
 * see  L. Fredman, R. Sedgewick, D. Sleator, and R. Tarjan,
 * The Pairing Heap: A New Form of Self-Adjusting Heap,
 * Algorithmica, 1986(1), pages 111-129.
 *
 * Converted from an Java implementation by Alexander Schwartz
 *
 * 07/30/2001
 * 07/01/2003 (updated)
 * 15/03/2004 (Copy comparision class)
 */


#include <functional>                 // for less<...>
#include <boost/concept_check.hpp>    // use the concept checking mechanism of boost
#include <iostream>


/**@class PairingHeap
 * @brief Implements a pairing heap
 *
 * @param Object  class of objects to store in heap
 * @param Cmp     class to compare objects
 *
 * The class @p Cmp should implement the 'Strict Weak Ordering' concept of the STL.
 * Hence, it also implements the 'Binary Predicate Concept' and therefore must
 * implement the following function:
 *
 * int operator() (const Object &, const Object &)
 *
 * 'std::less<Object>' is the default implementation, which can be used if the class Object
 * has implemented "<", i.e., @p Object implements the STL 'LessThanComparable' concept.
 */
template <typename Object, typename Cmp = std::less<Object> >
class PairingHeap
{
public:

  /**@class Item
   * @brief basic items in the heap for storing objects (of type @p Object)
   *
   * Each item keeps a copy of the object
   * (we take a copy since 'Object' is typically of small size).
   */
  class Item
  {
  public:
    friend class PairingHeap<Object, Cmp>;

    //! Constructor - gets a reference to obj, but creates a copy
    Item(Object &obj): parent(0), child(0), sibling(0), data(obj) {}

    //! Add the @p item as a child
    void addChild(Item *item)
    {
      item->sibling = child;
      item->parent  = this;
      this->child   = item;
    }

    //! Cut this item from its parent item
    void cutFromParent()
    {
      if (parent == 0)
	return;
      if (parent->child == this)
	parent->child = sibling;
      else
      {
	 Item *sister = parent->child;
	 while (sister->sibling != this)
	   sister = sister->sibling;
	 sister->sibling = sibling;
      }
      sibling = 0;
      parent = 0;
    }

  protected:
    Item *parent;
    Item *child;
    Item *sibling;
    Object data;
  };



 public:

  //! Construct an empty heap with default comparison class 'std::less<Object>'
  PairingHeap(bool check = false);

  /**@brief Constructor for an empty heap
   *
   * @param comp   an instance of the comparison function
   * @param check  whether a consistency check should be applied
   *
   * @p comp is now under the control of the heap and must NOT be deleted
   */
  PairingHeap(Cmp& comp, bool check = false);

  //! Destructor (delete comparison class instance)
  ~PairingHeap();

  //! returns the comparison class
  Cmp& getCmp() const { return comp_; }

  //! returns the number of elements in the heap
  int size() const { return length_; }

  //! checks whether the heap is empty
  bool isEmpty() const { return (root_ == 0); }

  //! returns the Object with minimal key
  Object getMin() const { return root_->data; }

  //! removes @p item from the heap
  void remove(Item* item);

  //! remove the item with minimum key
  Object removeMin();

  //! If the key of @p item has changed, this forces a new sorting.
  void decreaseKey(Item* item);

  //! add an object to the heap
  Item* add(Object &obj);

  //! Melt the current heap and @p P into a new heap
  void melt(PairingHeap<Object, Cmp>* P);


private:
  //! @name Internal functions
  //@{

  void delete_heap(Item *item);              //!< delete everything
  void checkInvariant();                     //!< check consistency
  void checkNode(Item *item, Item *parent);  //!< checks the nodes recursively
  void setRoot(Item *root);                  //!< sets a new root node
  Item* linkChildren(Item *parent);          //!< link the children below @p parent to the heap
  void pairChildren(Item *parent);           //!< pair the children below @p parent
  Item* link(Item *a, Item *b);              //!< link the given items
  //@}


  // ---------------------------------------------------------------------
  //! @name Variables
  //@{
  Cmp& comp_;
  Item* root_;
  int  length_;
  bool checkInvariant_;
  //@}
};




template <typename Object, typename Cmp>
PairingHeap<Object, Cmp>::PairingHeap(bool check) : comp_()
{
  root_ = 0;
  length_ = 0;
  checkInvariant_ = check;
}

template <typename Object, typename Cmp>
PairingHeap<Object, Cmp>::PairingHeap(Cmp& comp, bool check) : comp_(comp)
{
  root_ = 0;
  length_ = 0;
  checkInvariant_ = check;
}

template <typename Object, typename Cmp>
PairingHeap<Object, Cmp>::~PairingHeap()
{
  if (root_ != 0)
    delete_heap(root_);
}


// removes the item from the heap
template <typename Object, typename Cmp>
void PairingHeap<Object,Cmp>::remove(Item *item)
{
  if (item == root_)
  {
     removeMin();
     return;
  }
  item->cutFromParent();
  PairingHeap<Object,Cmp> *P = new PairingHeap<Object,Cmp>(comp_, checkInvariant_);
  P->length_ = 1;
  P->root_ = item;
  P->removeMin();
  melt(P);
}


// remove the item with minimum key
template <typename Object, typename Cmp>
Object PairingHeap<Object,Cmp>::removeMin()
{
  --length_;
  Object data = root_->data;
  Item *child = root_->child;
  if (child == 0)
  {
     delete root_;
     setRoot(0);
     if (checkInvariant_)
       checkInvariant();
     return data;
  }
  Item *save_root = root_;
  pairChildren(root_);
  setRoot(linkChildren(root_));
  delete save_root;
  if (checkInvariant_)
    checkInvariant();
  return data;
}


// if the key of the item has changed force a sorting
template <typename Object, typename Cmp>
void PairingHeap<Object,Cmp>::decreaseKey(Item* item)
{
  if (item == root_)
    return;

  item->cutFromParent();
  setRoot(link(root_, item));

  if (checkInvariant_)
    checkInvariant();
}


// add an Oject to the heap
template <typename Object, typename Cmp>
typename PairingHeap<Object,Cmp>::Item* PairingHeap<Object, Cmp>::add(Object &obj)
{
  ++length_;
  Item *item = new Item(obj);
  setRoot(link(root_, item));
  if (checkInvariant_)
    checkInvariant();
  return item;
}


// melt the current heap and Q into a new heap
template <typename Object, typename Cmp>
void PairingHeap<Object,Cmp>::melt(PairingHeap<Object,Cmp>* P)
{
  length_ += P->length_;
  setRoot(link(root_, P->root_));
  P->root_ = 0;
  if (checkInvariant_)
    checkInvariant();
}


// delete everything
template <typename Object, typename Cmp>
void PairingHeap<Object,Cmp>::delete_heap(Item* item)
{
  if (item->sibling != 0)
    delete_heap(item->sibling);
  if (item->child != 0)
    delete_heap(item->child);
  delete item;
}


// checks the representation invariant
template <typename Object, typename Cmp>
void PairingHeap<Object,Cmp>::checkInvariant()
{
  checkNode(root_, 0);
}


// checks the nodes recursively
template <typename Object, typename Cmp>
void PairingHeap<Object,Cmp>::checkNode(Item* item, Item* parent)
{
  if (item == 0)
    return;

  if (item->parent != parent)
  {
     std::cerr << "Error: node " << item << " has wrong parent (" << item->parent;
     std::cerr << " instead of " << parent << ")" << std::endl;
     exit(1);
  }
  if (item == parent)
  {
     std::cerr << "Error: node " << item << " is its own parent!" << std::endl;
     exit(1);
  }
  checkNode(item->sibling, parent);
  checkNode(item->child, item);
}


// sets the root
template <typename Object, typename Cmp>
void PairingHeap<Object,Cmp>::setRoot(Item* root)
{
  root_ = root;
  if (root_ != 0)
  {
     root_->sibling = 0;
     root_->parent  = 0;
  }
}


// link children
template <typename Object, typename Cmp>
typename PairingHeap<Object,Cmp>::Item *PairingHeap<Object,Cmp>::linkChildren(Item* parent)
{
  Item* child = parent->child;
  Item* myroot = child;

  if (child == 0)
    return myroot;

  child = child->sibling;
  myroot->sibling = 0;
  while (child != 0)
  {
     Item* next = child->sibling;
     child->sibling = 0;
     myroot = link(myroot, child);
     child = next;
  }
  return myroot;
}


// pair children
template <typename Object, typename Cmp>
void PairingHeap<Object,Cmp>::pairChildren(Item* parent)
{
  Item* firstChild = parent->child;
  if (firstChild == 0)
    return;

  parent->child = 0;

  Item* next = firstChild;
  do
  {
     Item* child1 = next;
     Item* child2 = next->sibling;

     if (child2 == 0)
     {
	parent->addChild(child1);
	return;
     }
     next = child2->sibling;
     child1->sibling = 0;
     child2->sibling = 0;

     parent->addChild(link(child1, child2));
  }
  while(next != 0);
}

// Link two items
template <typename Object, typename Cmp>
typename PairingHeap<Object,Cmp>::Item *PairingHeap<Object,Cmp>::link(Item* a, Item* b)
{
  boost::function_requires< boost::BinaryPredicateConcept<Cmp, Object, Object> >();

  if (a == 0)
    return b;
  if (b == 0)
    return a;

  // compare -> call "operator()"
  if ( comp_(a->data, b->data) )
  {
     a->addChild(b);
     return a;
  }

  b->addChild(a);
  return b;
}

#endif
