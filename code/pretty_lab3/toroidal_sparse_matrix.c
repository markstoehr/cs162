/** @toroidal_sparse_matrix.c
 *  @brief Function implementations for a toroidal sparse matrix using 
 *         doubly linked circularly linked lists
 *
 *  This contains the implementations for the sparse matrix over a torus
 *  represented using doubly linked circularly linked lists.
 *  We have implemented efficient algorithms for insertion
 *  and deletion which is special to the data structure.  Some of
 *  the algorithms are not optimally efficient, however,
 *  since as an intermediate data structure we use linked lists
 *  rather than priority queues or sorting 
 *  when iterating over elements (see matrix_entry_list.h). This
 *  is for clarity purposes in the program.
 *
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>. 
 *  Copyright (C) 2015
 *
 *  @author Mark Stoehr ()
 *  @date 2015-02-02
 *  @bug No known bugs.
 */


#include "toroidal_sparse_matrix.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


MatrixColumn *MatrixColumn_create(int index) {
  MatrixColumn *column = (MatrixColumn *) malloc(sizeof(MatrixColumn));
  if (!(column))
    return NULL;

  MatrixNode *head = MatrixNode_create(0);
  if (!(head)) {
    free(column);
    return NULL;
  }

  column->head = head;
  head->column = column;
  head->row = NULL;

  column->left = column;
  column->right = column;
  column->size = 0;
  column->index = index;
  column->matrix = NULL;

  return column;
}

MatrixRow *MatrixRow_create(int index) {
  MatrixRow *row = (MatrixRow *) malloc(sizeof(MatrixRow));
  if (!(row))
    return NULL;

  MatrixNode *head = MatrixNode_create(0);
  if (!(head)) {
    free(row);
    return NULL;
  }

  row->head = head;
  head->column = NULL;
  head->row = row;

  row->up = row;
  row->down = row;
  row->size = 0;
  row->index = index;
  row->matrix = NULL;

  return row;
}

void MatrixColumn_destroy(MatrixColumn *column) {
  if (!(column->head)) {
    free(column);
    return;
  }
  MatrixNode *node = (column->head)->down;
  while (node != column->head) {
    node = node->down;
    MatrixNode_destroy(node);
    
  }

  MatrixNode_destroy(column->head);
  
  // gracefully remove the column
  if (column->left &&
      ((column->left)->right == column))
    (column->left)->right = column->right;

  if (column->right &&
      ((column->left)->right == column))
    (column->right)->left = column->left;

  free(column);
  return;
  
}

void MatrixRow_destroy(MatrixRow *row) {
  if (!(row->head)) {
    free(row);
    return;
  }
  MatrixNode *node = (row->head)->right;
  while (node != row->head) {
    node = node->right;
    MatrixNode_destroy(node);
    
  }

  MatrixNode_destroy(row->head);
  
  // gracefully remove the row
  if (row->up &&
      ((row->up)->down == row))
    (row->up)->down = row->down;

  if (row->down &&
      ((row->up)->down == row))
    (row->down)->up = row->up;

  free(row);
  return;
  
}

MatrixNode *MatrixNode_create(int value) {
  MatrixNode *node = (MatrixNode *) malloc(sizeof(MatrixNode));
  if (!(node)) 
    return NULL;

  node->left = node;
  node->right = node;
  node->up = node;
  node->down = node;

  node->column = NULL;
  node->row = NULL;

  node->value = value;
  return node;

}

void MatrixNode_destroy(MatrixNode *node) {
  
  if (!(node))
    return;
  // if the left and right nodes are defined
  // and point to the node to be deleted fix the references
  if (node->left &&
      ((node->left)->right == node)) 
    (node->left)->right = node->right;
  
  if (node->right &&
      ((node->right)->left == node)) 
    (node->right)->left = node->left;

  if (node->up &&
      ((node->up)->down == node))
    (node->up)->down = node->down;

  if (node->down &&
      ((node->down)->up == node))
    (node->down)->up == node;

  if (node->column)
    --(node->column)->size;

  if (node->row)
    --(node->row)->size;

  if ((node->column)->matrix)
    --((node->column)->matrix)->nnz;

  free(node);
  return;
}

void MatrixNode_pop(MatrixNode *node) {
  
  if (!(node))
    return;
  // if the left and right nodes are defined
  // and point to the node to be deleted fix the references
  if (node->left &&
      ((node->left)->right == node)) 
    (node->left)->right = node->right;
  
  if (node->right &&
      ((node->right)->left == node)) 
    (node->right)->left = node->left;

  if (node->up &&
      ((node->up)->down == node))
    (node->up)->down = node->down;

  if (node->down &&
      ((node->down)->up == node))
    (node->down)->up == node;

  if (node->column)
    --(node->column)->size;

  if (node->row)
    --(node->row)->size;

  if ((node->column)->matrix)
    --((node->column)->matrix)->nnz;
  
  return;
}


MatrixNode *MatrixColumn_access(MatrixColumn *column,int value) {
  MatrixNode *node = (column->head)->down;
  while ( (node != column->head) &&
	  (node->value <= value)) // key invariant keep searching while less than
    node = node->down;

  return node->up;
}

MatrixNode *MatrixRow_access(MatrixRow *row,int value) {
  MatrixNode *node = (row->head)->right;
  while ( (node != row->head) &&
	  ((node->column)->index <= value)) // key invariant keep searching while less than
    node = node->right;

  return node->left;
}

MatrixNode *MatrixColumn_push_value(MatrixColumn *column,int value) {
  if (!(column))
    return NULL;

  Matrix *node = MatrixColumn_access(column,value);
  if (node->value == value)
    return node;

  //otherwise allocate a new node
  Matrix *new_node = MatrixNode_create(value);
  
  new_node->up = node;
  new_node->down = node->down;
  (node->down)->up = new_node;
  node->down = new_node;
  new_node->column = node->column;
  ++(node->column)->size;
  return new_node;
}

MatrixNode *MatrixRow_push_value(MatrixColumn *row,int value) {
  if (!(row))
    return NULL;

  Matrix *node = MatrixRow_access(row,value);
  if (node->value == value)
    return node;

  //otherwise allocate a new node
  Matrix *new_node = MatrixNode_create(value);
  
  new_node->left = node;
  new_node->right = node->right;
  (node->right)->left = new_node;
  node->right = new_node;
  new_node->row = node->row;
  ++(node->row)->size;
  return new_node;
}

Matrix *Matrix_create() {
  Matrix *matrix = (Matrix *)malloc(sizeof(matrix));

  if (!(matrix))
    return NULL;

  MatrixColumn *headcolumn = (MatrixColumn *)malloc(sizeof(MatrixColumn));
  if (!(headcolumn)) {
    free(matrix);
    return NULL;
  }

  MatrixRow *headrow = (MatrixRow *)malloc(sizeof(MatrixRow));
  if (!(headrow)) {
    free(matrix);
    free(headcolumn);
    return NULL;
  }

  headcolumn->matrix = matrix;
  headcolumn->index = 0;
  headcolumn->size = 0;
  headcolumn->head = NULL;
  headcolumn->left = headcolumn;
  headcolumn->right = headcolumn;

  headrow->matrix = matrix;
  headrow->index = 0;
  headrow->size = 0;
  headrow->head = NULL;
  headrow->up = headrow;
  headrow->down = headrow;

  

 
}

void Matrix_destroy(Matrix *matrix) {
  if (!(matrix))
    return;

  MatrixColumn *column = (matrix->headcolumn)->right;
  while (column != matrix->headcolumn) {
    column = column->right;
    MatrixColumn_destroy(column->left);
  }
  MatrixColumn_destroy(column);

  MatrixRow *row = (matrix->headrow)->down;
  while (row != matrix->headrow) {
    row = row->down;
    MatrixRow_destroy(row->up);
  }
  MatrixRow_destroy(row);

  free(matrix);

  return;
}

MatrixColumn *Matrix_access_column(Matrix *matrix,int index) {
  if (!(matrix))
    return NULL;

  MatrixColumn *column = (matrix->headcolumn)->right;

  while ( (column != matrix->headcolumn) &&
	  (column->index <= index))
    column = column->right;

  return column->left;
}

MatrixRow *Matrix_access_row(Matrix *matrix,int index) {
  if (!(matrix))
    return NULL;

  MatrixRow *row = (matrix->headrow)->down;
  while ( (row != matrix->headrow) &&
	  (row->index <= index))
    row = row->down;

  return row->up;
}

Matrix *Matrix_push_column(Matrix *matrix, MatrixColumn *column) {
  if (!(matrix)  ||  !(column))
    return NULL;

  if (column->matrix)
    return NULL;
  
  Matrix *prec_column = Matrix_access_column(matrix,column->index);
  // fail if you can't insert because a column is already there
  // with that index
  if (prec_column->index == column->index) 
    return NULL;

  // otherwise do simple insertion
  column->right = prec_column->right;
  column->left = prec_column;
  (prec_column->right)->left = column;
  prec_column->right = column;

  column->matrix = matrix;

  matrix->nnz += column->size;
  ++matrix->nnz_columns;

  // now we hook up to the rows
  MatrixNode *node = (column->head)->down;
  MatrixRow *row = (matrix->headrow)->down;
  MatrixNode *row_node;
  MatrixRow *new_row;
  /* @invariant row != matrix->headrow
   *            node != column->head
   * @invariant (row->up)->index <= node->value
   *            (node->up)->value <= row->index
   *
   * @invariant row->index > node->value
   *            implies there is no row
   *            r such that r->index == node->value
   * 
   *
   */
  while (  (row != matrix->headrow) &&
	   (node != column->head) ) {
    
    if (row->index == node->value) {
      row_node = MatrixRow_access(row,column->index);
      node->left = row_node;
      node->right = row_node->right;
      (row_node->right)->left = node;
      row_node->right = node;
    
      row = row->down;
      node = node->down;
    } else if ( row->index < node->value) 
      row = row->down;
    else  {// row->index > node->value
      // can derive from invariant that
      // `node` has no corresponding MatrixRow in `matrix`
      new_row = MatrixRow_create(node->value);
      (new_row->head)->right = node;
      (new_row->head)->left = node;
      node->right = new_row->head;
      node->left = new_row->head;

      new_row->matrix = matrix;
      ++new_row->size;
      
      ++matrix->nnz_rows;
      
      rowo
    }
    
    

  }
  
  return matrix;
}
