;REPORT
(define outlet-values
  (lambda (arg)
    (for-each
      (lambda (variable)
        (ti-menu-load-string
          (format #f "report surface-integrals area-weighted-avg 4 () ~a no quit"
          variable))
        (newline))
      '(totc_ni totc_mn totc_co totc_nh3 totc_na totc_so4 eqc_oh supersat ph smd m0 m1 m2 m3 l1 l2 w1 w2))))

;REGION. Create the regions
(define create-region
  (lambda (arg)
    (for-each
      (lambda (name in/out min-point max-point)
        (ti-menu-load-string
          (format #f "solve cell-registers add ~a type hexahedron inside? ~a min-point ~a max-point ~a quit"
          name in/out min-point max-point))
        (newline))
      '(region_metal region_nh3 region_oh region_0)
      '(yes yes no yes)
      '("0.002 0.001 0.041" "-0.002 -0.001 0.041" "-0.012 -0.002 0" "-0.002 -0.002 0")
      '("0.012 0.002 0.042" "-0.012 -0.002 0.042" "0.012 0.002 0.042" "0.002 0.002 0.042"))))

;PLANE. Create planes to see results.
(define create-plane-line
  (lambda (arg)
    (for-each
      (lambda (name method value)
        (ti-menu-load-string
          (format #f "surface plane-surface ~a ~a ~a quit" name method value))
        (newline))
      '(xy yz-x=0 yz-x=-0.0015 yz-x=0.0015 xz-y=0 xz-y=-0.0015 xz-y=0.0015)
      '(xy-plane yz-plane yz-plane yz-plane zx-plane zx-plane zx-plane)
      '(0.0415 0 -0.0015 0.0015 0 -0.0015 0.0015))

    ;LINE. Create lines to see results.
    (for-each
      (lambda (name start end)
        (ti-menu-load-string
          (format #f "surface line-surface ~a ~a ~a quit" name start end))
        (newline))
      '(line-outlet line-z=0.04125 line-z=0.0415 line-z=0.04175)
      '("0 0 0" "-0.002 0 0.04125" "-0.002 0 0.0415" "-0.002 0 0.04175")
      '("0 0 0.042" "0.002 0 0.04125" "0.002 0 0.0415" "0.002 0 0.04175"))))
      
;DISCRETIZATION SCHEME. Set the order of discretization scheme
(define schemeorder
  (lambda (order udsindex)
    (for-each
      (lambda (index)
        (ti-menu-load-string
          (format #f "solve set discretization-scheme uds-~a ~a" index order))
        (newline))
      udsindex)))

;UNDER-RELAXATION. Set the under-relaxation factors
(define relaxfactor
  (lambda (rf udsindex)
    (for-each
      (lambda (index)
        (ti-menu-load-string
          (format #f "solve set under-relaxation uds-~a ~a" index rf))
        (newline))
      udsindex)))

;SOURCE TERM. Set the source for the UDS. If we want to activate the source we have to write
;"1 no yes and the name of the source (conc or mom)", otherwise we have to write "0"
(define sourceterm
  (lambda (selection)
    (cond ( (equal? selection "none")
            (ti-menu-load-string
              (format #f
                (string-append
                  (do ((str "define boundary-conditions fluid fluid-00 no yes 0 0 0 0 0 0\n") (i 0 (+ i 1))) ((= i 13) str)
                      (define str (string-append str "0\n")))
                  "no no no 0. no 0. no 0. no 0 no 0 no 1 no no no no no\n"))))
          ( (equal? selection "species")
            (ti-menu-load-string
              (format #f
                (string-append
                  (do ((str "define boundary-conditions fluid fluid-00 no yes 0 0 0 0 0 0\n") (i 0 (+ i 1))) ((= i 6) str)
                      (define str (string-append str "1 no yes \"conc_source::libudf\"\n")))
                  (do ((str "") (i 0 (+ i 1))) ((= i 7) str)
                      (define str (string-append str "0\n")))
                  "no no no 0. no 0. no 0. no 0 no 0 no 1 no no no no no\n"))))
          ( (equal? selection "environments")
            (ti-menu-load-string
              (format #f
                (string-append
                  (do ((str "define boundary-conditions fluid fluid-00 no yes 0 0 0 0 0 0\n") (i 0 (+ i 1))) ((= i 6) str)
                      (define str (string-append str "0\n")))
                  (do ((str "") (i 0 (+ i 1))) ((= i 3) str)
                      (define str (string-append str "1 no yes \"env_source::libudf\"\n")))
                  (do ((str "") (i 0 (+ i 1))) ((= i 4) str)
                      (define str (string-append str "0\n")))
                  "no no no 0. no 0. no 0. no 0 no 0 no 1 no no no no no\n"))))
          ( (equal? selection "moments")
            (ti-menu-load-string
              (format #f
                (string-append
                  (do ((str "define boundary-conditions fluid fluid-00 no yes 0 0 0 0 0 0\n") (i 0 (+ i 1))) ((= i 9) str)
                      (define str (string-append str "0\n")))
                  (do ((str "") (i 0 (+ i 1))) ((= i 4) str)
                      (define str (string-append str "1 no yes \"mom_source::libudf\"\n")))
                  "no no no 0. no 0. no 0. no 0 no 0 no 1 no no no no no\n"))))
          ( (equal? selection "all")
            (ti-menu-load-string
              (format #f
                (string-append
                  (do ((str "define boundary-conditions fluid fluid-00 no yes 0 0 0 0 0 0\n") (i 0 (+ i 1))) ((= i 6) str)
                      (define str (string-append str "1 no yes \"conc_source::libudf\"\n")))
                  (do ((str "") (i 0 (+ i 1))) ((= i 3) str)
                      (define str (string-append str "1 no yes \"env_source::libudf\"\n")))
                  (do ((str "") (i 0 (+ i 1))) ((= i 4) str)
                      (define str (string-append str "1 no yes \"mom_source::libudf\"\n")))
                  "no no no 0. no 0. no 0. no 0 no 0 no 1 no no no no no\n")))))))

;PARAMETERS. Set the parameters of aggregation
  (define setparameters
    (lambda (names values)
      (for-each
        (lambda (name value)
          (if (not (rp-var-object name))
            (rp-var-define name value 'real #f)
            (rpsetvar name value))
          (newline))
        names
        values)))

(define initialize-precNMC
  (lambda (arg)
    ;CUSTOM FIELD FUNCTION. Define the custom field function to patch the udm-21
    (ti-menu-load-string
      (format #f "define custom-field-functions define \"diss-rate\" \"turb_diss_rate\" "))

    ;USER-DEFINED. Set UDSs, UDMs and UDF
    (for-each
      (lambda (type number options)
        (ti-menu-load-string
          (format #f "define user-defined ~a ~a ~a" type number options)))
      '(user-defined-scalars user-defined-memory "compiled-functions load" "function-hooks adjust")
      '(13 26 "" "")
      (list
        (do ((str "no no\n") (i 0 (+ i 1))) ((= i 13) str)
            (define str (string-append str "yes\n\"mass flow rate\"\n")))
        "\n" "\"libudf\"\n" "\n\"adjust::libudf\"\n"))

    ;DIFFUSIVITY. Set UDS diffusivity (0 for the moments)
    (ti-menu-load-string
      (format #f
        (string-append
          (do ((str "define materials change-create water-liquid water-liquid no no no no no no yes defined-per-uds\n") (i 0 (+ i 1))) ((= i 6) str)
              (define str (string-append str (format #f "~a user-defined \"conc_diffusivity::libudf\"\n" i))))
          (do ((str "") (i 6 (+ i 1))) ((= i 13) str)
              (define str (string-append str (format #f "~a user-defined \"env_diffusivity::libudf\"\n" i))))
          "-1 no \n")))

    ; SCHEME. Select SIMPLE algorithm
    (ti-menu-load-string
      (format #f "solve set p-v-coupling 20"))
    (newline)

    ;EQUATIONS. Set equations to solve
    (for-each
      (lambda (variable flag)
        (ti-menu-load-string
          (format #f "solve set equations ~a ~a" variable flag))
        (newline))
      '(flow kw)
      '(no no))

    (do ((i 0 (+ i 1))) ((> i 12))
      (ti-menu-load-string
        (format #f "solve set equations uds-~a yes" i))
        (newline))

    ;CONVERGENCE-CRITERIA. Set convergence criteria
    (newline)
    (ti-menu-load-string
      (format #f
        (do ((str "solve monitors residual convergence-criteria\n") (i 0 (+ i 1))) ((= i 13) str)
            (define str (string-append str "1e-20\n")))))

    ;Set scheme order
    (schemeorder 0 '(0 1 2 3 4 5 6 7 8 9 10 11 12))

    ;Set the under-relaxation factor
    (relaxfactor 1e-5 '(0 1 2 3 4 5 6 7 8 9 10 11 12))

    ;BOUNDARY
    ;Velocity-inlet.
    ;The numeric values are VELOCITY, GAUGE PRESSURE, TURBULENCE INTENSITY, HYDRAULIC DIAMETER,
    ;the concentrations and UDS 5, 6, 7, 8
    (for-each
      (lambda (inlet vel turb-intensity hd p1 p2 p3)
        (ti-menu-load-string
          (format #f "define boundary-conditions velocity-inlet ~a no no yes yes no ~a no 0. no no no yes ~a ~a no yes no yes no yes no yes no yes no yes no yes no yes no yes no yes no yes no yes no yes no 0. no 0. no 0. no 0. no 0. no 0. no ~a no ~a no ~a no 0. no 0. no 0. no 0. "
            inlet vel turb-intensity hd p1 p2 p3))
        (newline))
      '(inlet_metals inlet-nh3 inlet-oh01 inlet-oh02)
      (make-list 4 0.374)
      (make-list 4 7.727)
      (make-list 4 0.001)
      '(1 0 0 0)
      '(0 1 0 0)
      '(0 0 1 1))

    ;BOUNDARY
    ;Pressure-outlet. The numeric values are GAUGE PRESSURE, TURBULENCE INTENSITY, HYDRAULIC DIAMETER and the UDS 0-8
    (ti-menu-load-string
      (format #f "define boundary-conditions pressure-outlet outlet yes no 0. no yes no no no yes 1 0.002 yes yes yes yes yes yes yes yes yes yes yes yes yes no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 no 0 yes yes no no" ))

    ;Activate source terms
    (sourceterm "all")

    ;PATCH
    (do ((i 0 (+ i 1))) ((> i 12))
      (ti-menu-load-string
        (format #f "solve patch (fluid-00) () uds-~a no 0 " i))
      (newline))

    (do ((i 0 (+ i 1))) ((> i 24))
      (ti-menu-load-string
        (format #f "solve patch (fluid-00) () udm-~a no 0 " i))
      (newline))

    (for-each
      (lambda (region index)
        (ti-menu-load-string
          (format #f "solve patch () (~a) uds-~a no 1 " region index))
        (newline))
      '(region_metal region_nh3 region_oh)
      '(6 7 8))

    (ti-menu-load-string (format #f "solve patch (fluid-00) () udm-25 yes diss-rate "))
    (newline)
    
    (ti-menu-load-string (format #f "solve patch () (region_0) udm-24 no 1.0 "))))
