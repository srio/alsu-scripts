; (c) Kenneth A. Goldberg  http://goldberg.lbl.gov  KAGoldberg(at)lbl.gov
; Lawrence Berkeley National Laboratory
; 09/05/2011, based on linear_gsfit1.pro
; IDL files: linear_2dgsfit1.pro
;
;  u are the values we want to fit (the target array)
;    on the same domain as b[i].a
;
; b = linear_2dpolynomials1(N, Ny=Ny, Nz=12, /double)  ;--- create the basis. Ny is the optional domain width
; v = linear_2dgsfit1(wavefront, b, Nz=8)              ;--- fit a wavefront to 8+1 polynomials (order 8)
; wf_fit = linear_basis1(v, b)                         ;--- show the fit wavefront phase
;
; 09/05/2011 created this to handle rectangular apertures, but general enough that
; it should work with any orthogonal polynomials.

; u is an array that you want to fit to the basis, b
; b is the structure containing the basis polynomials.
; NZ can be inferred from the size of b, or specified by the user. It's the # of fitting terms.
; mask is a point-by-point weighting with the same size as u
;
; IMPORTANT: Never switch from using a mask to not using a mask with the same basis functions.
;    Orthogonalization is performed USING the mask if it's present, or not. So be consistent.
; IMPORTANT: input array u, b[i].a basis functions, and mask must all be of the same array size
;
; Note, the scalar values b[i].total_squared are calculated in the orthonormalization.
;   When the polynomials are orthonormal, total_squared is 1.0 by definition.
;   When they are merely orthogonal, they take a positive-definite value.
;   They are calculated with  b[i].total_squared = total( mask * b[i].a^2, /double ), or
;     b[i].total_squared = total( b[i].a^2, /double )  depending on whether or not a mask (weighting)
;     is being used.
;

function linear_2dgsfit1, u, b, Nz=Nz, mask=mask

    ;--- Either use a user-specified Nz, or infer it from the number of basis functions in the structure, b
    ;    Make sure we don't exceed the number of basis functions.
  Nz = (defined(Nz) ? Nz : n_elements(b) ) < n_elements(b)

    ;--- This will be the vector of fit coefficients. Use Nz+1 output terms by convnetion. 
    ;    Calculations will be double or single precision based on the type of b[i].a
    ;
  v = (size(b[0].a, /type) EQ 5) ? dblarr(Nz) : fltarr(Nz)   ;--- output vector of coefficients

  if defined(mask) then begin   ;--- Use the mask weighting function.

    for i=0,Nz-1 do $           ;    Project onto the basis, normalize properly
      v[i] = total(mask * u * b[i].a) / b[i].total_squared  ;--- here, total_squared respects the mask

  endif else begin              ;--- Use no mask weighting function.

    for i=0,Nz-1 do $           ;    Project onto the basis, normalize properly
      v[i] = total(u * b[i].a) / b[i].total_squared

  endelse

return, v  ;--- return the output coefficients in a vector of size Nz+1
end
