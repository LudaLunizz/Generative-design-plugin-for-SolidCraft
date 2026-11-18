//! Which GPU adapters wgpu sees, and whether they offer 64-bit floats in shaders.
fn main() {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor { backends: wgpu::Backends::VULKAN, ..wgpu::InstanceDescriptor::new_without_display_handle() });
    for a in pollster::block_on(instance.enumerate_adapters(wgpu::Backends::VULKAN)) {
        let i = a.get_info();
        let f = a.features();
        let l = a.limits();
        println!("{} | {:?} | driver {} {} | f64 {} | max storage buffer {} MiB | max workgroup invocations {}", i.name, i.device_type, i.driver, i.driver_info, f.contains(wgpu::Features::SHADER_F64), l.max_storage_buffer_binding_size / (1 << 20), l.max_compute_invocations_per_workgroup);
    }
}
