(define pic_interval 0.00833333)
(define current (rpgetvar 'flow-time))

(define (run-images)
	(if (< (+ pic_interval current) (rpgetvar 'flow-time)) 
		(begin 
			(set! current (+ pic_interval current))
			(ti-menu-load-string "dis set-w 2")
			(ti-menu-load-string "dis con solids vof 0 0.6")
			(ti-menu-load-string "dis h-c Verify_CA-0.884mm-solids_vf-%f.jpeg")
			(ti-menu-load-string "dis set-w 3")
			(ti-menu-load-string "dis con mix log_temp_dif -4 4")
			(ti-menu-load-string "dis h-c Verify_CA-0.884mm-log_temp_dif-%f.jpeg")
			(ti-menu-load-string "dis set-w 4")
			(ti-menu-load-string "dis con mix conversion 0 7")
			(ti-menu-load-string "dis h-c Verify_CA-0.884mm-conversion-%f.jpeg")
			(ti-menu-load-string "dis set-w 5")
			(ti-menu-load-string "dis con air uds-0-scalar 0 1")
			(ti-menu-load-string "dis h-c Verify_CA-0.884mm-air_scalar-%f.jpeg")
			(ti-menu-load-string "dis set-w 6")
			(ti-menu-load-string "dis con solids uds-1-scalar 0 1")
			(ti-menu-load-string "dis h-c Verify_CA-0.884mm-solids_scalar-%f.jpeg")
		)
	)
	(if (> (rpgetvar 'flow-time) (+ pic_interval current)) 
		(begin 
			(set! current (rpgetvar 'flow-time))
		)
	)
	(if (< (rpgetvar 'flow-time) current) 
		(begin 
			(set! current (rpgetvar 'flow-time))
		)
	)
)
