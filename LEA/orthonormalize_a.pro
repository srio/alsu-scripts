; (c) Kenneth A. Goldberg  http://goldberg.lbl.gov  KAGoldberg(at)lbl.gov
; Lawrence Berkeley National Laboratory
; 1/3/20
; IDL files: orthonormalize_a.pro
;
;--- This function was written to orthonormalize the APS AXO data influence functions on 2020-01-03,
;    but it is general and can orthonormalize any-dimensional basis functions with or without weighting
;
;    The input argument, a, is a structure variable defined like this example
;       a = replicate({a:dblarr(500), total_squared:0d}, 38)
;    In this way, we access a[0].a, a[1].a, a[2].a as basis functions.
;    In the original intention of this routine, a[i].a will contain the influence functions
;      from the APS AXO mirror (measured in October 2019). The basis functions of a[i].a are
;      not necessarily orthogonal or orthonormal.
;    To facilitate fitting, we add a constant a[i].a = 1, and a linear term a[i].a = x,
;      to the basis.
;
;   The optional 'mask' keyword is a user-supplied weighting function.
;   It should have the same dimensions as a[i].a .
;   If mask is not given, we assume uniform weighting across the domain.
;
;   The 'matrix' keyword returns a NxN transformation matrix from the input basis, a, 
;     to the output basis, b.
;
;   2020-01-27, I changed the normalization so that
;   total(mask*b[i].a^2) = 1.0 for all i
;
 
function orthonormalize_a, a, mask=mask, matrix=matrix

    ;--- a has to be defined like this
    ;    a = replicate({a:x, total_squared:0d}, Nz+1)

  b = a  ;--- copy so it's the same size as the original. b becomes the output, the orthonormal basis.
  Nz = n_elements(a)-1
  matrix = dblarr(Nz+1, Nz+1)  ;--- this will be returned in the keyword

  if isa(mask) then begin  ;--- do this code if there is a mask

      ;--- This is a Gram-Schmidt process
      ;    One by one, each new elements is projected onto the existing basis and made orthogonal.
    for i=1,Nz do begin
      for j=0,i-1 do begin
        matrix[i,j] = (total(mask * b[i].a * b[j].a, /double) / total( mask * b[j].a^2, /double ))
        b[i].a -= b[j].a * matrix[i,j]
      endfor
      b[i].a /= sqrt(total(mask * b[i].a^2, /double))  ;--- normalize each element to rms = 1.0
    endfor

      ;--- compute and store the total_squared for each. We use this in future normalization to save time
    for i=0,Nz do $
      b[i].total_squared = total( mask * b[i].a^2, /double )  ;--- must be 1.0

  endif else begin   ;--- do this code below if there is no mask
  
      ;--- This is a Gram-Schmidt process
      ;    One by one, each new elements is projected onto the existing basis and made orthogonal.      
    for i=1,Nz do begin
      for j=0,i-1 do begin
        matrix[i,j] = (total(b[i].a * b[j].a, /double) / total( b[j].a^2, /double ))
        b[i].a -= b[j].a * matrix[i,j]
      endfor
      b[i].a /= sqrt(total(b[i].a^2, /double ))  ;--- normalize each element to rms = 1.0
    endfor

      ;--- compute and store the total_squared for each. We use this in future normalization to save time
    for i=0,Nz do $
      b[i].total_squared = total( b[i].a^2, /double )  ;--- must be 1.0
  endelse

return, b  ;--- return the structure containing the basis.
end