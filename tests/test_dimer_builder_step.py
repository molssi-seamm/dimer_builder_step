#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `dimer_builder_step` package."""

import pytest  # noqa: F401
import dimer_builder_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = dimer_builder_step.DimerBuilder()
    assert (
        str(type(result)) == "<class 'dimer_builder_step.dimer_builder.DimerBuilder'>"
    )
