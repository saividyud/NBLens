"""
GPU node stats script for Unity cluster testing.
Outputs VRAM, node name, device info, and verifies CuPy can access the allocated GPU.
Run via: sbatch "Unity/Bash Scripts/GPU Tests/gpu_test.sh"
"""

import platform
import sys


def main():
    print("=" * 60)
    print("GPU Node Statistics (Unity Cluster Test)")
    print("=" * 60)

    # Host / node info (available without GPU)
    print("\n--- Host Info ---")
    print(f"Node name:        {platform.node()}")
    print(f"Machine:          {platform.machine()}")
    print(f"System:           {platform.system()}")
    print(f"Processor:        {platform.processor()}")

    # GPU info via CuPy
    print("\n--- GPU Info ---")
    try:
        import cupy as cp

        num_devices = cp.cuda.runtime.getDeviceCount()
        print(f"Visible GPUs:     {num_devices}")

        if num_devices == 0:
            print("WARNING: No GPUs detected. Check SLURM --gres=gpu:1 and CUDA_VISIBLE_DEVICES.")
            sys.exit(1)

        for i in range(num_devices):
            device = cp.cuda.Device(i)
            free_mem, total_mem = device.mem_info

            # Device name from CUDA runtime (cudaDeviceProp.name)
            try:
                props = cp.cuda.runtime.getDeviceProperties(i)
                name_val = getattr(props, "name", None)
                if name_val is not None:
                    device_name = (name_val.split(b"\x00")[0].decode() if isinstance(name_val, bytes) else str(name_val)).strip()
                else:
                    device_name = "Unknown"
            except Exception:
                device_name = "Unknown"

            print(f"\n  GPU {i}:")
            print(f"    Name:           {device_name}")
            print(f"    Total VRAM:     {total_mem / 1024**3:.2f} GB")
            print(f"    Free VRAM:      {free_mem / 1024**3:.2f} GB")
            print(f"    Used VRAM:      {(total_mem - free_mem) / 1024**3:.2f} GB")
            print(f"    Compute cap:    {device.compute_capability}")
            print(f"    PCI Bus ID:     {device.pci_bus_id}")

        # Quick CuPy sanity check
        print("\n--- CuPy Sanity Check ---")
        x = cp.array([1.0, 2.0, 3.0])
        y = cp.sum(x)
        print(f"  CuPy test (sum of [1,2,3] on GPU): {float(y)}")
        print("  CuPy GPU compute: OK")

    except ImportError as e:
        print(f"ERROR: CuPy not installed or not in environment: {e}")
        print("  Activate your .venv and ensure cupy-cuda12x (or matching CUDA) is installed.")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print("\n" + "=" * 60)
    print("GPU test completed successfully.")
    print("=" * 60)


if __name__ == "__main__":
    main()
