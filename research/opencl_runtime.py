"""Shared OpenCL device/context/buffer infrastructure, extracted from the MIT Chua project.

No model kernel is compiled here. Model subclasses compile their own paired kernels.
"""
import numpy as np
import pyopencl as cl

class OpenCLRuntime:
    def __init__(self, device_filter="", precision="float64"):
        devices = [d for p in cl.get_platforms() for d in p.get_devices()
                   if d.type & cl.device_type.GPU and device_filter.lower() in d.name.lower()]
        if not devices:
            raise RuntimeError("No matching OpenCL GPU. Install a GPU driver with OpenCL support. "
                               "This program does not silently fall back to CPU.")
        self.device = devices[0]
        if precision not in ("float32", "float64"):
            raise ValueError("precision must be float32 or float64")
        if precision == "float64" and not self.device.double_fp_config:
            raise RuntimeError("This GPU does not support float64; explicitly select float32.")
        self.dtype = np.dtype(precision)
        self.scalar = self.dtype.type
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context, properties=cl.command_queue_properties.PROFILING_ENABLE)
        self.info = {"name": self.device.name.strip(), "vendor": self.device.vendor,
                     "driver": self.device.driver_version, "opencl": self.device.version,
                     "precision": precision, "global_memory_bytes": self.device.global_mem_size,
                     "backend": "PyOpenCL", "fast_math": False}

    def buffer(self, array, readonly=False):
        flags = cl.mem_flags.READ_ONLY if readonly else cl.mem_flags.READ_WRITE
        return cl.Buffer(self.context, flags | cl.mem_flags.COPY_HOST_PTR, hostbuf=array)

