JULIA_DEPOT_PATH := $(shell pwd)/.julenv

define julia_env
	julia --project=$(JULIA_DEPOT_PATH) $(1)
endef

juliaenv:
	@$(call julia_env,$(ARGS))

instantiateenv: 
	@$(call julia_env,-e 'using Pkg; Pkg.instantiate()')

.PHONY: juliaenv instantiateenv 
