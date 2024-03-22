include .env
#NOTE: prefiex commands with icn will work only if there's a .env containing 
#the appropiate credentials 

JULIA_DEPOT_PATH := $(shell pwd)/.julenv

define julia_env
	julia --project=$(JULIA_DEPOT_PATH) $(1)
endef

define icn_julia_env
	$(ICN_JULIA_BIN) --project=$(ICN_JULIA_DEPOT_PATH) $(1)
endef

juliaenv:
	@$(call julia_env,$(ARGS))

instantiatenv: 
	@$(call julia_env,-e 'using Pkg; Pkg.instantiate()')

icnjuliaenv:
	@$(call icn_julia_env,$(ARGS))

icninstantiatenv: 
	@$(call icn_julia_env,-e 'using Pkg; Pkg.instantiate()')

.PHONY: juliaenv instantiatenv icnjuliaenv icninstantiatenv 

