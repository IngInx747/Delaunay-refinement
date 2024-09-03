#ifndef OPENMESH_PROPERTY_WRAPPER_HH
#define OPENMESH_PROPERTY_WRAPPER_HH

#include <sstream>
#include <stdexcept>
#include <string>
#include <OpenMesh/Core/Utils/HandleToPropHandle.hh>
#include <OpenMesh/Core/Mesh/ArrayKernel.hh>

namespace OpenMesh {

namespace Property2 {

template <typename PROPTYPE, typename MeshT = int>
class PropertyManager
{
public:

    using Self            = PropertyManager<PROPTYPE, MeshT>;
    using Value           = typename PROPTYPE::Value;
    using value_type      = typename PROPTYPE::value_type;
    using Handle          = typename PROPTYPE::Handle;
    using Reference       = typename PROPTYPE::reference;
    using ConstReference  = typename PROPTYPE::const_reference;

protected:

    template <typename PropertyManager2, typename PropHandleT>
    struct StorageT;

    // definition for other Mesh Properties
    template <typename PropertyManager2, typename PropHandleT>
    struct StorageT
    {
        static void initialize(PropertyManager<PROPTYPE, MeshT>& pm, const Value& initial_value )
        { pm.set_range(pm.mesh_.template all_elements<Handle>(), initial_value); }

        static void copy(const PropertyManager& from, PropertyManager2& to)
        { from.copy_to(from.mesh_.template all_elements<Handle>(), to, to.mesh_.template all_elements<Handle>()); }

        static void swap(PropertyManager& lhs, PropertyManager2& rhs)
        {
            std::swap(lhs.mesh_.property(lhs.prop_).data_vector(), rhs.mesh_.property(rhs.prop_).data_vector());
            // resize the property to the correct size
            lhs.mesh_.property(lhs.prop_).resize(lhs.mesh_.template n_elements<Handle>());
            rhs.mesh_.property(rhs.prop_).resize(rhs.mesh_.template n_elements<Handle>());
        }

        static ConstReference access_property_const(ArrayKernel& mesh, const PROPTYPE& prop_handle, const Handle& handle)
        { return mesh.property(prop_handle, handle); }

        static Reference access_property(ArrayKernel& mesh, const PROPTYPE& prop_handle, const Handle& handle)
        { return mesh.property(prop_handle, handle); }
    };

    // specialization for Mesh Properties
    template <typename PropertyManager2>
    struct StorageT<PropertyManager2, MPropHandleT<Value>>
    {
        static void initialize(PropertyManager<PROPTYPE, MeshT>& pm, const Value& initial_value )
        { pm() = initial_value; }

        static void copy(const PropertyManager<PROPTYPE, MeshT>& from, PropertyManager2& to)
        { *to = *from; }

        static void swap(PropertyManager<PROPTYPE, MeshT>& from, PropertyManager2& to)
        { std::swap(*to, *from); }

        static ConstReference access_property_const(ArrayKernel& mesh, const PROPTYPE& prop_handle, const Handle&)
        { return mesh.property(prop_handle); }

        static Reference access_property(ArrayKernel& mesh, const PROPTYPE& prop_handle, const Handle&)
        { return mesh.property(prop_handle); }
    };

    using Storage = StorageT<Self, PROPTYPE>;

public:

    /// Constructor
    PropertyManager() = delete;

    /// Constructor
    /// Asks for a property with name propname and creates one if none exists. Lifetime is not managed.
    /// @param mesh The mesh on which to create the property.
    /// @param propname The name of the property.
    PropertyManager(ArrayKernel& _mesh, const char *propname) : mesh_(_mesh), retain_(true), name_(propname)
    {
        if (!mesh_.get_property_handle(prop_, propname))
            mesh_ref().add_property(prop_, propname);
    }

    /// Constructor
    /// Asks for a property with name propname and creates one if none exists. Lifetime is not managed.
    /// @param initial_value If the proeprty is newly created, it will be initialized with initial_value.
    ///        If the property already existed, nothing is changes.
    /// @param mesh The mesh on which to create the property.
    /// @param propname The name of the property.
    PropertyManager(const Value& initial_value, ArrayKernel& _mesh, const char *propname) : mesh_(_mesh), retain_(true), name_(propname)
    {
        if (!mesh_.get_property_handle(prop_, propname))
        {
            mesh_ref().add_property(prop_, propname);
            Storage::initialize(*this, initial_value);
        }
    }

    /// Constructor
    /// Create an anonymous property. Lifetime is managed.
    /// @param mesh The mesh on which to create the property.
    explicit PropertyManager(const ArrayKernel& _mesh) : mesh_(_mesh), retain_(false), name_("")
    {
        mesh_ref().add_property(prop_, name_);
    }

    /// Constructor
    /// Create an anonymous property. Lifetime is managed.
    /// @param initial_value The property will be initialized with initial_value.
    /// @param mesh The mesh on which to create the property.
    PropertyManager(const Value& initial_value, const ArrayKernel& _mesh) : mesh_(_mesh), retain_(false), name_("")
    {
        mesh_ref().add_property(prop_, name_);
        Storage::initialize(*this, initial_value);
    }

    /// Constructor
    /// Create a wrapper around an existing property. Lifetime is not managed.
    /// @param mesh The mesh on which to create the property.
    /// @param property_handle Handle to an existing property that should be wrapped.
    PropertyManager(ArrayKernel& mesh, PROPTYPE property_handle) : mesh_(mesh), prop_(property_handle), retain_(true), name_()
    {}

    /// Constructor
    /// Create a wrapper around an existing property. Lifetime is not managed.
    /// @param mesh The mesh on which to view the property.
    /// @param property_handle Handle to an existing property that should be wrapped.
    PropertyManager(const ArrayKernel& mesh, PROPTYPE property_handle) : mesh_(mesh), prop_(property_handle), retain_(true), name_()
    {}

    PropertyManager(const PropertyManager& rhs) : mesh_(rhs.mesh_), prop_(), retain_(rhs.retain_), name_(rhs.name_)
    {
        if (rhs.retain_) // named property -> create a property manager referring to the same
        {
            prop_ = rhs.prop_;
        }
        else // unnamed property -> create a property manager refering to a new property and copy the contents
        {
            mesh_ref().add_property(prop_, name_);
            Storage::copy(rhs, *this);
        }
    }

    /// Create property manager referring to a copy of the current property.
    /// This can be used to explicitely create a copy of a named property. The cloned property
    /// will be unnamed.
    PropertyManager clone()
    {
        PropertyManager result(mesh_);
        Storage::copy(*this, result);
        return result;
    }

    PropertyManager& operator=(const PropertyManager& rhs)
    {
        if (&mesh_ == &rhs.mesh_ && prop_ == rhs.prop_)
            ; // nothing to do
        else
            Storage::copy(rhs, *this);
        return *this;
    }

    /// Move constructor. Transfers ownership (delete responsibility).
    PropertyManager(PropertyManager &&rhs) : mesh_(rhs.mesh_), prop_(rhs.prop_), retain_(rhs.retain_), name_(rhs.name_)
    {
        if (!rhs.retain_) // only invalidate unnamed properties
            rhs.prop_.invalidate();
    }

    /// Move assignment. Transfers ownership (delete responsibility).
    PropertyManager& operator=(PropertyManager&& rhs)
    {
        if ((&mesh_ != &rhs.mesh_) || (prop_ != rhs.prop_))
        {
            if (rhs.retain_)
            {
                // retained properties cannot be invalidated. Copy instead
                Storage::copy(rhs, *this);
            }
            else
            {
                // swap the data stored in the properties
                Storage::swap(rhs, *this);
                // remove the property from rhs
                rhs.mesh_ref().remove_property(rhs.prop_);
                // invalidate prop_
                rhs.prop_.invalidate();
            }
        }
        return *this;
    }

    /// Create a property manager for the supplied property and mesh.
    /// If the property doesn't exist, it is created. In any case,
    /// lifecycle management is disabled.
    /// @see makePropertyManagerFromExistingOrNew
    static PropertyManager createIfNotExists(ArrayKernel &mesh, const char *propname)
    { return PropertyManager(mesh, propname); }

    /// Like createIfNotExists() with two parameters except, if the property
    /// doesn't exist, it is initialized with the supplied value over
    /// the supplied range after creation. If the property already exists,
    /// this method has the exact same effect as the two parameter version.
    /// Lifecycle management is disabled in any case.
    /// @see makePropertyManagerFromExistingOrNew
    template<typename PROP_VALUE, typename ITERATOR_TYPE>
    static PropertyManager createIfNotExists(
        ArrayKernel &mesh,
        const char *propname,
        const ITERATOR_TYPE &begin,
        const ITERATOR_TYPE &end,
        const PROP_VALUE &init_value)
    {
        const bool exists = propertyExists(mesh, propname);
        PropertyManager pm(mesh, propname, exists);
        pm.retain();
        if (!exists) pm.set_range(begin, end, init_value);
        return std::move(pm);
    }

    /// Like createIfNotExists() with two parameters except, if the property
    /// doesn't exist, it is initialized with the supplied value over
    /// the supplied range after creation. If the property already exists,
    /// this method has the exact same effect as the two parameter version.
    /// Lifecycle management is disabled in any case.
    /// @see makePropertyManagerFromExistingOrNew
    template<typename PROP_VALUE, typename ITERATOR_RANGE>
    static PropertyManager createIfNotExists(
        ArrayKernel &mesh,
        const char *propname,
        const ITERATOR_RANGE &range,
        const PROP_VALUE &init_value)
    {
        return createIfNotExists(mesh, propname, range.begin(), range.end(), init_value);
    }


    /// Access the value of the encapsulated mesh property.
    /// Example:
    /// @code
    /// GraphMesh m;
    /// auto description = getOrMakeProperty<void, std::string>(m, "description");
    /// *description = "This is a very nice mesh.";
    /// @endcode
    /// @note This method is only used for mesh properties.
    typename PROPTYPE::reference& operator*()
    {
        return mesh_ref().mproperty(prop_)[0];
    }

    /// Access the value of the encapsulated mesh property.
    /// Example:
    /// @code
    /// GraphMesh m;
    /// auto description = getProperty<void, std::string>(m, "description");
    /// std::cout << *description << std::endl;
    /// @endcode
    /// @note This method is only used for mesh properties.
    typename PROPTYPE::const_reference& operator*() const
    {
        return mesh_ref().mproperty(prop_)[0];
    }

    /// Enables convenient access to the encapsulated property.
    /// For a usage example see this class' documentation.
    /// @param handle A handle of the appropriate handle type. (I.e. \p VertexHandle for \p VPropHandleT, etc.)
    inline typename PROPTYPE::reference operator[] (Handle handle)
    {
        return mesh_ref().property(prop_, handle);
    }

    /// Enables convenient access to the encapsulated property.
    /// For a usage example see this class' documentation.
    /// @param handle A handle of the appropriate handle type. (I.e. \p VertexHandle for \p VPropHandleT, etc.)
    inline typename PROPTYPE::const_reference operator[] (const Handle& handle) const
    {
        return mesh_ref().property(prop_, handle);
    }

    /// Enables convenient access to the encapsulated property.
    /// For a usage example see this class' documentation.
    /// @param handle A handle of the appropriate handle type. (I.e. \p VertexHandle for \p VPropHandleT, etc.)
    inline typename PROPTYPE::reference operator() (const Handle& handle = Handle())
    {
        // return mesh_ref().property(prop_, handle);
        return Storage::access_property(mesh_ref(), prop_, handle);
    }

    /// Enables convenient access to the encapsulated property.
    /// For a usage example see this class' documentation.
    /// @param handle A handle of the appropriate handle type. (I.e. \p VertexHandle for \p VPropHandleT, etc.)
    inline typename PROPTYPE::const_reference operator() (const Handle& handle = Handle()) const
    {
        // return mesh_ref().property(prop_, handle);
        return Storage::access_property_const(mesh_ref(), prop_, handle);
    }

    /// 
    /// Conveniently set the property for an entire range of values.
    /// 
    /// Examples:
    /// \code
    /// MeshT mesh;
    /// PropertyManager<VPropHandleT<double>> distance(
    ///     mesh, "distance.plugin-example.i8.informatik.rwth-aachen.de");
    /// distance.set_range(
    ///     mesh.vertices_begin(), mesh.vertices_end(),
    ///     std::numeric_limits<double>::infinity());
    /// \endcode
    /// or
    /// \code
    /// MeshT::VertexHandle vh;
    /// distance.set_range(
    ///     mesh.vv_begin(vh), mesh.vv_end(vh),
    ///     std::numeric_limits<double>::infinity());
    /// \endcode
    /// 
    /// @param begin Start iterator. Needs to dereference to HandleType.
    /// @param end End iterator. (Exclusive.)
    /// @param value The value the range will be set to.
    /// 
    template<typename HandleTypeIterator, typename PROP_VALUE>
    void set_range(HandleTypeIterator begin, HandleTypeIterator end, const PROP_VALUE &value)
    {
        for (; begin != end; ++begin)
            (*this)[*begin] = value;
    }

    /// 
    /// Conveniently transfer the values managed by one property manager
    /// onto the values managed by a different property manager.
    /// 
    /// @param begin Start iterator. Needs to dereference to HandleType. Will
    /// be used with this property manager.
    /// @param end End iterator. (Exclusive.) Will be used with this property
    /// manager.
    /// @param dst_propmanager The destination property manager.
    /// @param dst_begin Start iterator. Needs to dereference to the
    /// HandleType of dst_propmanager. Will be used with dst_propmanager.
    /// @param dst_end End iterator. (Exclusive.)
    /// Will be used with dst_propmanager. Used to double check the bounds.
    /// 
    template<typename HandleTypeIterator, typename PropertyManager2, typename HandleTypeIterator2>
    void copy_to(
        HandleTypeIterator begin,
        HandleTypeIterator end,
        PropertyManager2 &dst_propmanager,
        HandleTypeIterator2 dst_begin,
        HandleTypeIterator2 dst_end) const
    {
        for (; begin != end && dst_begin != dst_end; ++begin, ++dst_begin)
            dst_propmanager[*dst_begin] = (*this)[*begin];
    }

    template<typename RangeType, typename PropertyManager2, typename RangeType2>
    void copy_to(
        const RangeType &range,
        PropertyManager2 &dst_propmanager,
        const RangeType2 &dst_range) const
    {
        copy_to(range.begin(), range.end(), dst_propmanager, dst_range.begin(), dst_range.end());
    }

    /// 
    /// Copy the values of a property from a source range to
    /// a target range. The source range must not be smaller than the
    /// target range.
    /// 
    /// @param prop_name Name of the property to copy. Must exist on the
    /// source mesh. Will be created on the target mesh if it doesn't exist.
    /// 
    /// @param src_mesh Source mesh from which to copy.
    /// @param src_range Source range which to copy. Must not be smaller than
    /// dst_range.
    /// @param dst_mesh Destination mesh on which to copy.
    /// @param dst_range Destination range.
    /// 
    template<typename RangeType, typename RangeType2>
    static void copy(
        const char *prop_name,
        ArrayKernel &src_mesh, const RangeType  &src_range,
        ArrayKernel &dst_mesh, const RangeType2 &dst_range)
    {
        typedef PropertyManager<PROPTYPE> DstPM;
        DstPM dst(DstPM::createIfNotExists(dst_mesh, prop_name));
        typedef PropertyManager<PROPTYPE> SrcPM;
        SrcPM src(src_mesh, prop_name, true);
        src.copy_to(src_range, dst, dst_range);
    }

    /// Mark whether this property should be stored when mesh is written to a file
    /// @param _persistence Property will be stored iff _persistence is true
    void set_persistent(bool _persistence = true)
    {
        mesh_.property(getRawProperty()).set_persistent(_persistence);
    }

    ~PropertyManager()
    {
        deleteProperty();
    }

    // swap the data stored in the properties
    inline void swap(PropertyManager &rhs)
    {
        Storage::swap(rhs, *this);
    }

    static bool propertyExists(const ArrayKernel &mesh, const char *propname)
    {
        PROPTYPE dummy;
        return mesh.get_property_handle(dummy, propname);
    }

    inline bool isValid() const { return prop_.is_valid(); }

    inline operator bool() const { return isValid(); }

    inline const PROPTYPE &getRawProperty() const { return prop_; }

    inline const std::string &getName() const { return name_; }

    inline MeshT &getMesh() const { return dynamic_cast<const MeshT&>(mesh_); }

    template <typename MeshType>
    inline MeshType &getMesh() const { return dynamic_cast<const MeshType&>(mesh_); }

protected:

    void deleteProperty()
    {
        if (!retain_ && prop_.is_valid())
        {
            mesh_ref().remove_property(prop_);
        }
    }

    // Property changing is a non-constant operation
    // Something needs to be re-designed in kernel
    inline ArrayKernel &mesh_ref() const
    {
        return const_cast<ArrayKernel&>(mesh_);
    }

protected:

    const ArrayKernel &mesh_;

    PROPTYPE prop_;

    bool retain_;

    std::string name_;
};

#if 0

template <typename PropertyT>
class ConstPropertyViewer
{
public:

    using Value             = typename PropertyT::Value;
    using value_type        = typename PropertyT::value_type;
    using Handle            = typename PropertyT::Handle;
    using ConstReference    = typename PropertyT::const_reference;

protected:

    template <typename PropHandleT>
    struct StorageT;

    // definition for other Mesh Properties
    template <typename PropHandleT>
    struct StorageT
    {
        static ConstReference access_property_const(const ArrayKernel& mesh, const PropertyT& prop_handle, const Handle& handle)
        { return mesh.property(prop_handle, handle); }
    };

    // specialization for Mesh Properties
    template <>
    struct StorageT<MPropHandleT<Value>>
    {
        static ConstReference access_property_const(const ArrayKernel& mesh, const PropertyT& prop_handle, const Handle&)
        { return mesh.property(prop_handle); }
    };

    using Storage = StorageT<PropertyT>;

public:

    ConstPropertyViewer(const ArrayKernel& mesh, const PropertyT& property_handle)
    : mesh(mesh), prop(property_handle) {}

    inline const typename PropertyT::const_reference operator[] (const Handle& handle)
    { return Storage::access_property_const(mesh, prop, handle); }

    inline const typename PropertyT::const_reference operator() (const Handle& handle = Handle())
    { return Storage::access_property_const(mesh, prop, handle); }

protected:

    const ArrayKernel &mesh;

    PropertyT prop;
};

template<typename ElementT, typename T>
inline Property2::ConstPropertyViewer<typename HandleToPropHandle<ElementT, T>::type> getProperty(const ArrayKernel &mesh, const char *propname)
{
    typename HandleToPropHandle<ElementT, T>::type prop;

    if (!mesh.get_property_handle(prop, propname))
    {
        std::ostringstream oss;
        oss << "Requested property handle \"" << propname << "\" does not exist.";
        throw std::runtime_error(oss.str());
    }

    return Property2::ConstPropertyViewer<decltype(prop)>(mesh, prop);
}

#endif

} // namespace Property2

/// @relates PropertyManager
/// 
/// Creates a new property whose lifetime is limited to the current scope.
/// 
/// Used for temporary properties. Shadows any existing properties of
/// matching name and type.
/// 
/// Example:
/// @code
/// GraphMesh m;
/// {
///     auto is_quad = makeTemporaryProperty<FaceHandle, bool>(m);
///     for (auto& fh : m.faces()) {
///         is_quad[fh] = (m.valence(fh) == 4);
///     }
///     // The property is automatically removed from the mesh at the end of the scope.
/// }
/// @endcode
/// 
/// @param mesh The mesh on which the property is created
/// @tparam ElementT Element type of the created property, e.g. VertexHandle, HalfedgeHandle, etc.
/// @tparam T Value type of the created property, e.g., \p double, \p int, etc.
/// @returns A PropertyManager handling the lifecycle of the property
/// 
template<typename ElementT, typename T>
inline Property2::PropertyManager<typename HandleToPropHandle<ElementT, T>::type> makeTemporaryProperty(ArrayKernel &mesh)
{
    return Property2::PropertyManager<typename HandleToPropHandle<ElementT, T>::type>(mesh);
}

/// @relates PropertyManager
/// 
/// Tests whether a property with the given element type, value type, and name is
/// present on the given mesh.
/// 
/// * Example:
/// @code
/// GraphMesh m;
/// if (hasProperty<FaceHandle, bool>(m, "is_quad")) {
///     // We now know the property exists: getProperty won't throw.
///     auto is_quad = getProperty<FaceHandle, bool>(m, "is_quad");
///     // Use is_quad here.
/// }
/// @endcode
/// 
/// @param mesh The mesh in question
/// @param propname The property name of the expected property
/// @tparam ElementT Element type of the expected property, e.g. VertexHandle, HalfedgeHandle, etc.
/// @tparam T Value type of the expected property, e.g., \p double, \p int, etc.
/// @tparam MeshT Type of the mesh. Can often be inferred from \p mesh
/// 
template<typename ElementT, typename T>
inline bool hasProperty(const ArrayKernel &mesh, const char *propname)
{
    typename HandleToPropHandle<ElementT, T>::type ph;
    return mesh.get_property_handle(ph, propname);
}

/// @relates PropertyManager
/// 
/// Obtains a handle to a named property.
/// 
/// Example:
/// @code
/// GraphMesh m;
/// {
///     try {
///         auto is_quad = getProperty<FaceHandle, bool>(m, "is_quad");
///         // Use is_quad here.
///     }
///     catch (const std::runtime_error& e) {
///         // There is no is_quad face property on the mesh.
///     }
/// }
/// @endcode
/// 
/// @pre Property with the name \p propname of matching type exists.
/// @throws std::runtime_error if no property with the name \p propname of
/// matching type exists.
/// @param mesh The mesh on which the property is created
/// @param propname The name of the created property
/// @tparam ElementT Element type of the created property, e.g. VertexHandle, HalfedgeHandle, etc.
/// @tparam T Value type of the created property, e.g., \p double, \p int, etc.
/// @returns A PropertyManager wrapping the property
/// 
template<typename ElementT, typename T>
inline Property2::PropertyManager<typename HandleToPropHandle<ElementT, T>::type> getProperty(ArrayKernel &mesh, const char *propname)
{
    if (!hasProperty<ElementT, T>(mesh, propname))
    {
        std::ostringstream oss;
        oss << "Requested property handle \"" << propname << "\" does not exist.";
        throw std::runtime_error(oss.str());
    }
    return Property2::PropertyManager<typename HandleToPropHandle<ElementT, T>::type>(mesh, propname);
}

/// @relates PropertyManager
/// 
/// Obtains a handle to a named property if it exists or creates a new one otherwise.
/// 
/// Used for creating or accessing permanent properties.
/// 
/// Example:
/// @code
/// GraphMesh m;
/// {
///     auto is_quad = getOrMakeProperty<FaceHandle, bool>(m, "is_quad");
///     for (auto& fh : m.faces()) {
///         is_quad[fh] = (m.valence(fh) == 4);
///     }
///     // The property remains on the mesh after the end of the scope.
/// }
/// {
///     // Retrieve the property from the previous scope.
///     auto is_quad = getOrMakeProperty<FaceHandle, bool>(m, "is_quad");
///     // Use is_quad here.
/// }
/// @endcode
/// 
/// @param mesh The mesh on which the property is created
/// @param propname The name of the created property
/// @tparam ElementT Element type of the created property, e.g. VertexHandle, HalfedgeHandle, etc.
/// @tparam T Value type of the created property, e.g., \p double, \p int, etc.
/// @returns A PropertyManager wrapping the property
/// 
template<typename ElementT, typename T>
inline Property2::PropertyManager<typename HandleToPropHandle<ElementT, T>::type> getOrMakeProperty(ArrayKernel &mesh, const char *propname)
{
    return Property2::PropertyManager<typename HandleToPropHandle<ElementT, T>::type>::createIfNotExists(mesh, propname);
}

/// @relates PropertyManager
/// @pre Property with the name \p propname of matching type exists.
/// @throws std::runtime_error if no property with the name \p propname of matching type exists.
/// @param mesh The mesh on which the property is created
/// @param propname The name of the created property
/// @tparam ElementT Element type of the created property, e.g. VertexHandle, HalfedgeHandle, etc.
/// @tparam T Value type of the created property, e.g., \p double, \p int, etc.
/// @returns A ConstPropertyViewer wrapping the property
/// 
template<typename ElementT, typename T>
inline Property2::PropertyManager<typename HandleToPropHandle<ElementT, T>::type> getProperty(const ArrayKernel &mesh, const char *propname)
{
    typename HandleToPropHandle<ElementT, T>::type prop;

    if (!mesh.get_property_handle(prop, propname))
    {
        std::ostringstream oss;
        oss << "Requested property handle \"" << propname << "\" does not exist.";
        throw std::runtime_error(oss.str());
    }

    return Property2::PropertyManager<typename HandleToPropHandle<ElementT, T>::type>(mesh, prop);
}

template<typename HandleT, typename ValueT>
inline const ValueT *getPropertyPtr(const ArrayKernel &mesh, const char *propname)
{
    typename HandleToPropHandle<HandleT, ValueT>::type prop;
    if (!mesh.get_property_handle(prop, propname)) return nullptr;
    return mesh.property(prop).data();
}

template<typename HandleT, typename ValueT>
inline void removeProperty(ArrayKernel &mesh, const char *propname)
{
    typename HandleToPropHandle<HandleT, ValueT>::type prop;
    if (mesh.get_property_handle(prop, propname)) mesh.remove_property(prop);
}

} // namespace OpenMesh

#endif