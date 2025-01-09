(define pic_interval 0.0002)
(define current (rpgetvar 'flow-time))

(define (export-images)
	(if (< (+ pic_interval current) (rpgetvar 'flow-time)) 
		(begin 
			(set! current (+ pic_interval current))
			;; xy-view
			(ti-menu-load-string "/views/restore-view xy")
			(ti-menu-load-string "/display/set/contours/surfaces xsec_yx 0 ,")
			(ti-menu-load-string "/display/contour solid vof 0 0.63")
			(ti-menu-load-string "/display/save-picture xy_solid_vof-%5.4f.png")
			(ti-menu-load-string "/display/contour fluid velocity 0 0.5")
			(ti-menu-load-string "/display/save-picture xy_fluid_vel-%5.4f.png")
			(ti-menu-load-string "/display/contour solid velocity 0 0.5")
			(ti-menu-load-string "/display/save-picture xy_solid_vel-%5.4f.png")
			;; yz-view
			(ti-menu-load-string "/views/restore-view yz")
			(ti-menu-load-string "/display/set/contours/surfaces xsec_yz 0 ,")
			(ti-menu-load-string "/display/contour solid vof 0 0.63")
			(ti-menu-load-string "/display/save-picture yz_solid_vof-%5.4f.png")
			(ti-menu-load-string "/display/contour fluid velocity 0 0.5")
			(ti-menu-load-string "/display/save-picture yz_fluid_vel-%5.4f.png")
			(ti-menu-load-string "/display/contour solid velocity 0 0.5")
			(ti-menu-load-string "/display/save-picture yz_solid_vel-%5.4f.png")
			;; xyz-view
			(ti-menu-load-string "/views/restore-view xyz")
			(ti-menu-load-string "/display/set/contours/surfaces xsec_yx xsec_yz ,")
			(ti-menu-load-string "/display/contour solid vof 0 0.63")
			(ti-menu-load-string "/display/save-picture xyz_solid_vof-%5.4f.png")
			(ti-menu-load-string "/display/contour fluid velocity 0 0.5")
			(ti-menu-load-string "/display/save-picture xyz_fluid_vel-%5.4f.png")
			(ti-menu-load-string "/display/contour solid velocity 0 0.5")
			(ti-menu-load-string "/display/save-picture xyz_solid_vel-%5.4f.png")
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



