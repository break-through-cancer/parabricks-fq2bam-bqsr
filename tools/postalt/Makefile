IMAGE ?= ghcr.io/break-through-cancer/parabricks-fq2bam-bqsr-postalt:0.2.1
PLATFORM ?= linux/amd64

.PHONY: build test

build:
	docker build --platform $(PLATFORM) -t $(IMAGE) .

test: build
	docker run --rm --platform $(PLATFORM) -v $(CURDIR)/tests:/tests:ro $(IMAGE) bash /tests/test_postalt_apply.sh
