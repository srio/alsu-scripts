; (c) Kenneth A. Goldberg  http://goldberg.lbl.gov  KAGoldberg(at)lbl.gov
; Lawrence Berkeley National Laboratory
; 8/03/09
; IDL files: linear_basis1.pro

; When the linear polynomial basis functions, b[i].a, are already calculated,
; using linear_polynomials1.pro, or a similar routine, we can calculate an output vector
; from an input set of coefficients, v
;
; b = linear_polynomials1(N, Nz=12, /double)  ;--- create the basis
; v = linear_gsfit1(wavefront, b, Nz=8)       ;--- fit a wavefront to 8 polynomials
; wf_fit = linear_basis1(v, b)                ;--- show the fit wavefront phase
;
; v is a vector of coefficients
; b is a structure containing the basis functions
; This routine does not worry about the size of the arrays

function linear_basis1, v, b

    ;--- Create an output array filled with zeros to start
    ;    The array will be single or double precision depending on the basis function type.
    ;
  y = b[0].a * 0
  
    ;--- If the vector is too long, complain, but don't fail
  if n_elements(v) GT n_elements(b) $
    then message, 'Mismatch between the number of coefficients and number of polynomials', /INFORMATIONAL
  
    ;--- choose an output number of elements N that does not exceed the number of basis functions
  N = n_elements(v) < n_elements(b)
  
    ;--- Add up the vectors. We use a for loop because the basis functions are contained
    ;    in a structure making it hard to create a giant matric multiplication.
    ;
  for i=0,N-1 do $
    if (v[i] NE 0.) then y += v[i] * b[i].a   ;--- we skip zeros to save time.

return, y  ;--- the output array will have the same size as b[i].a elements
end